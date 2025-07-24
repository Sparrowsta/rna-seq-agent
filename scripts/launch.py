import json
import subprocess
import os
import sys
import re
from collections import defaultdict
import argparse
import concurrent.futures
from datetime import datetime


def setup_directories():
    """
    创建必要的文件夹用以存储过程中产生的文件
    """
    data_dir = '/data'
    
    sub_dirs = [
        "srr_files",
        "fastq_files",
        "fastp",
        "bam",
        "featurecounts",
        "genomes",
        "work",
        ".nxf_home"
    ]

    try:
        os.makedirs(data_dir, exist_ok=True)
        print(f"Directory created or already exists: {data_dir}")

        for sub_dir in sub_dirs:
            full_path = os.path.join(data_dir, sub_dir)
            os.makedirs(full_path, exist_ok=True)
            print(f"Directory created or already exists: {full_path}")
            
        return True
    except OSError as e:
        print(f"Error creating directories: {e}")
        return False


def download_and_unzip(url, dest_folder):
    """Downloads a file from a URL, saves it, and unzips if it's a .gz file."""
    
    file_name = os.path.basename(url.split('?')[0])
    gz_path = os.path.join(dest_folder, file_name)
    unzipped_path = gz_path.replace(".gz", "")

    if os.path.exists(unzipped_path):
        print(f"File already exists, skipping download: {unzipped_path}")
        return unzipped_path

    print(f"Downloading {url} to {gz_path}...")
    try:
        subprocess.run(["wget", "-O", gz_path, url], check=True, capture_output=True, text=True)
        
        print(f"Unzipping {gz_path}...")
        subprocess.run(["gunzip", "-f", gz_path], check=True)
        
        print(f"Successfully prepared: {unzipped_path}")
        return unzipped_path
    except subprocess.CalledProcessError as e:
        print(f"Error during download/unzip for {url}.", file=sys.stderr)
        print(f"Stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'wget' or 'gunzip' command not found. Please ensure they are installed and in your PATH.", file=sys.stderr)
        sys.exit(1)

def process_single_srr(srr_id, fastq_dir, seq_mode):
    """Processes a single SRR ID: download, convert, and compress."""
    try:
        print(f"Processing SRR ID: {srr_id}")
        
        final_fastq_files = [os.path.join(fastq_dir, f"{srr_id}.fastq.gz")]
        if not seq_mode:
            final_fastq_files = [
                os.path.join(fastq_dir, f"{srr_id}_1.fastq.gz"),
                os.path.join(fastq_dir, f"{srr_id}_2.fastq.gz")
            ]

        if all(os.path.exists(f) for f in final_fastq_files):
            return f"{srr_id}: SKIPPED - FASTQ.gz files already exist."

        sra_dir = os.path.join(os.path.dirname(fastq_dir), "srr_files")
        os.makedirs(sra_dir, exist_ok=True)
        sra_path = os.path.join(sra_dir, srr_id)
        srr_file = os.path.join(sra_path, f"{srr_id}.sra")

        if not os.path.exists(srr_file):
            print(f"Downloading {srr_id}...")
            subprocess.run(
                ["bash", "-c", f"source activate sra_env && prefetch -O {sra_dir} {srr_id} -f ALL -H 1"],
                check=True, text=True, capture_output=True, timeout=1800
            )
        else:
            print(f"SRA file for {srr_id} already exists.")

        print(f"Converting {srr_id} to FASTQ and compressing...")
        
        if seq_mode:
            command = f"""
            source activate sra_env && \\
            fasterq-dump --stdout --split-spot --skip-technical -t {sra_dir} {srr_file} | gzip > {final_fastq_files[0]}
            """
            subprocess.run(["bash", "-c", command], check=True, text=True, capture_output=True, timeout=1800)
        else:
            dump_command = f"source activate sra_env && fasterq-dump --split-files -t {sra_dir} -O {fastq_dir} {srr_file}"
            subprocess.run(["bash", "-c", dump_command], check=True, text=True, capture_output=True, timeout=1800)
            
            gzip_r1_command = f"gzip {os.path.join(fastq_dir, srr_id + '_1.fastq')}"
            gzip_r2_command = f"gzip {os.path.join(fastq_dir, srr_id + '_2.fastq')}"
            subprocess.run(["bash", "-c", gzip_r1_command], check=True, text=True, capture_output=True)
            subprocess.run(["bash", "-c", gzip_r2_command], check=True, text=True, capture_output=True)

        return f"{srr_id}: SUCCESS"

    except subprocess.CalledProcessError as e:
        return f"{srr_id}: FAILED - {e.stderr}"
    except subprocess.TimeoutExpired:
        return f"{srr_id}: FAILED - Timeout"
    except Exception as e:
        return f"{srr_id}: FAILED - An unexpected error occurred: {e}"

def prepare_fastq_files(srr_list_path, fastq_dir, seq_mode, max_workers=4):
    """Reads SRR list, processes them in parallel, and verifies all outputs in parallel."""
    print(f"\n--- Preparing FASTQ Files (in parallel with max {max_workers} workers) ---")
    if not os.path.isfile(srr_list_path):
        print(f"Error: SRR list file not found at {srr_list_path}", file=sys.stderr)
        return False

    with open(srr_list_path, 'r') as f:
        srr_ids = [line.strip() for line in f if line.strip()]

    if not srr_ids:
        print(f"Error: SRR list file is empty: {srr_list_path}", file=sys.stderr)
        return False

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_srr = {executor.submit(process_single_srr, srr_id, fastq_dir, seq_mode): srr_id for srr_id in srr_ids}
        for future in concurrent.futures.as_completed(future_to_srr):
            srr_id = future_to_srr[future]
            try:
                result = future.result()
                print(f"Result for {srr_id}: {result}")
            except Exception as exc:
                print(f'{srr_id} generated an exception: {exc}')

    print("--- FASTQ File Preparation Complete ---")
    return True

def prepare_pipeline_command(params: dict) -> list:
    """
    Prepares all files and constructs the Nextflow command list.
    This function is intended to be called by the MCP Server.
    It does NOT execute the command.
    """
    print("--- Preparing Pipeline Command with Provided Parameters ---")
    print(json.dumps(params, indent=4, default=str))

    # 1. Setup directories
    if not setup_directories():
        raise RuntimeError("Failed to set up directories.")

    # 2. Prepare FASTQ files if SRR list is provided
    fastq_dir = params.get("fastq_dir") or '/data/fastq_files'
    os.makedirs(fastq_dir, exist_ok=True)
    
    if params.get("srr_list"):
        if not prepare_fastq_files(params["srr_list"], fastq_dir, params["seq_mode"]):
            raise RuntimeError("Failed to prepare FASTQ files.")
    
    # 3. Prepare genome files
    glob_pattern = '*.fastq.gz' if params["seq_mode"] else '*_{1,2}.fastq.gz'
    params["reads"] = os.path.join(fastq_dir, glob_pattern)

    genome_data = parse_genome_file('/data/genome.txt')
    if not genome_data:
        raise RuntimeError("Could not parse genome data.")
        
    selected_species = params["species"]
    selected_version = params["genome_version"]
    
    version_data = genome_data.get(selected_species, {}).get(selected_version, {})
    source = version_data.get('source', 'remote')

    if source == 'local':
        print(f"Using local genome for {selected_species} {selected_version}")
        params["fasta"] = version_data.get("fasta")
        params["gtf"] = version_data.get("gtf")
    else:
        fasta_url = version_data.get("fasta")
        gtf_url = version_data.get("gtf")
        genome_base_dir = os.path.join('/data', 'genomes', selected_species, selected_version)
        os.makedirs(genome_base_dir, exist_ok=True)
        if fasta_url:
            params["fasta"] = download_and_unzip(fasta_url, genome_base_dir)
        if gtf_url:
            params["gtf"] = download_and_unzip(gtf_url, genome_base_dir)
    
    # 4. Construct Nextflow command
    log_dir = '/data/.nextflow/log'
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'nextflow.log')

    nextflow_params = []
    for key, value in params.items():
        if isinstance(value, bool):
            if value:
                nextflow_params.append(f"--{key}")
        elif value is not None:
            nextflow_params.extend([f"--{key}", str(value)])

    command = [
        "nextflow",
        "-log", log_file,
        "run", "/app/nextflow/main.nf",
        "-profile", "inside_container,standard",
        "-resume",
        "-work-dir", "/data/work"
    ] + nextflow_params

    print("\nPrepared command:")
    print(" ".join(command))
    
    return command
