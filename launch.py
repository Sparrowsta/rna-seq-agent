import json
import subprocess
import os
import sys
import re
from collections import defaultdict
import argparse
import concurrent.futures

def get_args():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="NGS Pipeline Launcher")
    parser.add_argument("--srr_list", type=str, help="Path to the SRR list file (for non-interactive mode).")
    parser.add_argument("--seq_mode", action="store_true", help="Flag for single-end sequencing.")
    parser.add_argument("--species", type=str, help="Species to process (for non-interactive mode),etc. human,mouse,other.")
    parser.add_argument("--genome_version", type=str, help="Genome version, e.g., hg38 (for non-interactive mode).")
    parser.add_argument("--fastq_dir", type=str, help="Path to the directory with local fastq files.")
    return parser.parse_args()

def setup_directories(base_dir):
    """
    Creates necessary directories for the NGS pipeline within the 'data' folder.
    Returns True on success, False on failure.
    """
    data_dir = os.path.join(base_dir, 'data')
    
    sub_dirs = [
        "srr_files",
        "fastq_files",
        "fastp",
        "bam",
        "featurecounts",
        "genomes"
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

def parse_genome_file(filepath, base_dir):
    """
    Parses the genome.txt file and scans local directories to return a structured dictionary.
    Categorizes genomes into 'human', 'mouse', and 'other'.
    Handles both nested (species/version) and flat (species/version.fa) local genome structures.
    """
    genome_data = defaultdict(lambda: defaultdict(dict))

    # 1. Parse genome.txt for downloadable genomes
    pattern = re.compile(r'([\w.-]+)\s+(fasta|gtf)\s*:\s*(https?://\S+)')
    try:
        with open(filepath, 'r') as f:
            for line in f:
                match = pattern.match(line.strip())
                if match:
                    filename, file_type, url = match.groups()
                    # Use the basename of the file as the version identifier
                    version = re.sub(r'\.(fa|fasta|gtf)$', '', filename)
                    
                    if 'hg' in version:
                        species = 'human'
                    elif 'mm' in version:
                        species = 'mouse'
                    else:
                        species = 'other'
                    
                    genome_data[species][version][file_type] = url
                    genome_data[species][version]['source'] = 'remote'

    except FileNotFoundError:
        print(f"Warning: Genome file not found at {filepath}. Only local genomes will be available.")

    # 2. Scan local genomes directory
    genomes_dir = os.path.join(base_dir, 'data', 'genomes')
    if os.path.exists(genomes_dir):
        for species in os.listdir(genomes_dir):  # e.g., human, mouse, other
            species_dir = os.path.join(genomes_dir, species)
            if not os.path.isdir(species_dir):
                continue

            # Accumulate files found directly within the species directory (e.g., other/danio.fa)
            flat_files = defaultdict(dict)

            for item_name in os.listdir(species_dir):
                item_path = os.path.join(species_dir, item_name)
                
                # Case 1: Item is a subdirectory (e.g., human/hg38/, other/danio_rerio/)
                if os.path.isdir(item_path):
                    version = item_name
                    fasta_file = next((f for f in os.listdir(item_path) if f.endswith(('.fa', '.fasta'))), None)
                    gtf_file = next((f for f in os.listdir(item_path) if f.endswith('.gtf')), None)
                    
                    if fasta_file and gtf_file and version not in genome_data[species]:
                        genome_data[species][version]['fasta'] = os.path.join(item_path, fasta_file)
                        genome_data[species][version]['gtf'] = os.path.join(item_path, gtf_file)
                        genome_data[species][version]['source'] = 'local'
                
                # Case 2: Item is a file (relevant for 'other' species, e.g., other/danio.fa)
                elif os.path.isfile(item_path) and species == 'other':
                    version = item_name.split('.')[0]
                    if item_name.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
                        flat_files[version]['fasta'] = item_path
                    elif item_name.endswith(('.gtf', '.gtf.gz')):
                        flat_files[version]['gtf'] = item_path

            # Process the collected flat files for the 'other' species
            if species == 'other':
                for version, files in flat_files.items():
                    if 'fasta' in files and 'gtf' in files and version not in genome_data[species]:
                        genome_data[species][version]['fasta'] = files['fasta']
                        genome_data[species][version]['gtf'] = files['gtf']
                        genome_data[species][version]['source'] = 'local'

    return json.loads(json.dumps(genome_data))

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
            fasterq-dump --stdout --split-spot --skip-technical {srr_file} | gzip > {final_fastq_files[0]}
            """
            subprocess.run(["bash", "-c", command], check=True, text=True, capture_output=True, timeout=1800)
        else:
            dump_command = f"source activate sra_env && fasterq-dump --split-files -O {fastq_dir} {srr_file}"
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

def verify_fastq_file(f_path):
    """Verifies a single FASTQ file for existence, size, and integrity."""
    if not os.path.exists(f_path) or os.path.getsize(f_path) < 100:
        print(f"ERROR: Final file missing or too small: {f_path}")
        return False
    try:
        subprocess.run(["gzip", "-t", f_path], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        print(f"ERROR: Final file is corrupt: {f_path}")
        return False

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

    print("\n--- Final Verification of All FASTQ Files ---")
    all_expected_files = []
    for srr_id in srr_ids:
        if seq_mode:
            all_expected_files.append(os.path.join(fastq_dir, f"{srr_id}.fastq.gz"))
        else:
            all_expected_files.append(os.path.join(fastq_dir, f"{srr_id}_1.fastq.gz"))
            all_expected_files.append(os.path.join(fastq_dir, f"{srr_id}_2.fastq.gz"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(verify_fastq_file, all_expected_files)
        if not all(results):
            print("\nERROR: One or more FASTQ files are missing or corrupt. Please check the logs.")
            print("Pipeline will not be started.")
            return False

    print("--- FASTQ File Preparation and Verification Complete ---")
    return True

def handle_non_interactive_mode(args):
    """Handles non-interactive mode."""
    print("\n--- Running in Non-Interactive Mode ---")
    errors = []
    if args.srr_list and args.fastq_dir:
        errors.append("Error: --srr_list and --fastq_dir cannot be used at the same time.")
    elif not args.srr_list and not args.fastq_dir:
        errors.append("Error: Either --srr_list or --fastq_dir must be provided.")
    
    if not args.species:
        errors.append("Error: --species is required.")
    if not args.genome_version:
        errors.append("Error: --genome_version is required.")

    if errors:
        print("\nParameter validation failed:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    params = {
        "srr_list": args.srr_list,
        "fastq_dir": args.fastq_dir,
        "seq_mode": args.seq_mode,
        "species": args.species,
        "genome_version": args.genome_version,
        "seq_type": "rna-seq"
    }
    
    print("Parameters loaded successfully:")
    if params['srr_list']:
        print(f"  SRR list file: {params['srr_list']}")
    if params['fastq_dir']:
        print(f"  FASTQ directory: {params['fastq_dir']}")
    print(f"  Single-end: {params['seq_mode']}")
    print(f"  Species: {params['species']}")
    print(f"  Genome Version: {params['genome_version']}")
    
    return params, False

def get_user_choice(prompt, options):
    """Gets user choice from a list of options."""
    while True:
        choice = input(prompt)
        if choice in options:
            return choice
        else:
            print("Invalid selection. Please try again.")

def get_data_source(current_dir):
    """Gets the data source from the user."""
    print("\nSelect the source of your FASTQ files:")
    print("  1: Download and process from an SRR list file (e.g., SRR_list.txt)")
    print("  2: Use existing FASTQ files in the 'fastq_files' directory")
    
    choice = get_user_choice("Enter number (1-2): ", ["1", "2"])
    
    if choice == "1":
        srr_list_path = input("Enter the path to your SRR list file (default: SRR_list.txt): ") or "SRR_list.txt"
        if not os.path.isfile(srr_list_path):
            print(f"Error: File not found at {srr_list_path}", file=sys.stderr)
            sys.exit(1)
        
        with open(srr_list_path, 'r') as f:
            if not any(line.strip() for line in f):
                print(f"Error: SRR list file is empty: {srr_list_path}", file=sys.stderr)
                sys.exit(1)
                
        return {"srr_list": srr_list_path}, True
    else:
        fastq_dir = os.path.join(current_dir, 'data', 'fastq_files')
        if not os.path.exists(fastq_dir) or not os.listdir(fastq_dir):
            print(f"Error: The directory '{fastq_dir}' is empty or does not exist.", file=sys.stderr)
            print("Please add your FASTQ files to this directory or choose to download from SRR.", file=sys.stderr)
            sys.exit(1)
        print("Using local FASTQ files. The SRR download step will be skipped.")
        return {"srr_list": None}, False

def get_seq_mode():
    """Gets the sequencing mode from the user."""
    print("\nIs the sequencing data single-end or paired-end?")
    print("  1: Single-End")
    print("  2: Paired-End")
    choice = get_user_choice("Enter number (1-2): ", ["1", "2"])
    return choice == "1"

def get_genome_selection(genome_data):
    """Gets the genome selection from the user."""
    all_species = list(genome_data.keys())
    species_options = {str(i+1): s for i, s in enumerate(all_species)}

    print("\nSelect the species category:")
    for key, value in species_options.items():
        print(f"  {key}: {value.capitalize()}")
    
    species_choice = get_user_choice(f"Enter number (1-{len(species_options)}): ", species_options.keys())
    selected_species = species_options[species_choice]

    version_options = {str(i+1): v for i, v in enumerate(genome_data[selected_species].keys())}
    print(f"\nSelect genome version for {selected_species.capitalize()}:")
    for key, value in version_options.items():
        source = genome_data[selected_species][value].get('source', 'remote')
        print(f"  {key}: {value} ({source})")
        
    version_choice = get_user_choice(f"Enter number (1-{len(version_options)}): ", version_options.keys())
    selected_version = version_options[version_choice]
    
    return selected_species, selected_version

def handle_interactive_mode(current_dir):
    """Handles interactive mode."""
    print("\n--- Running in Interactive Mode ---")
    params = {"seq_type": "rna-seq"}
    
    source_params, download_srr = get_data_source(current_dir)
    params.update(source_params)
    
    params["seq_mode"] = get_seq_mode()
    
    genome_data = parse_genome_file('genome.txt', current_dir)
    if not genome_data:
        sys.exit(1)
        
    species, genome_version = get_genome_selection(genome_data)
    params["species"] = species
    params["genome_version"] = genome_version
            
    return params, download_srr

def confirm_and_run(params, download_srr, is_interactive, current_dir):
    """Confirm parameters and run the pipeline."""
    if is_interactive:
        print("\n--- Configuration Summary ---")
        data_source_str = 'Local FASTQ files' if not download_srr else f"SRR List ({params.get('srr_list')})"
        print(f"  Data Source: {data_source_str}")
        print(f"  Sequencing: {'Single-End' if params.get('seq_mode') else 'Paired-End'}")
        print(f"  Species: {params.get('species')}")
        print(f"  Genome Version: {params.get('genome_version')}")
        print("-----------------------------")
        
        confirm = (input("Proceed with this configuration? [Y/n]: ").strip() or "y").lower()
        if confirm != 'y':
            print("Operation cancelled by user.")
            sys.exit(0)

    fastq_dir = params.get("fastq_dir") or os.path.join(current_dir, 'data', 'fastq_files')
    os.makedirs(fastq_dir, exist_ok=True)
    
    if params.get("srr_list"):
        if not prepare_fastq_files(params["srr_list"], fastq_dir, params["seq_mode"]):
            sys.exit(1)
    
    relative_fastq_dir = os.path.relpath(fastq_dir, current_dir)
    glob_pattern = '*.fastq.gz' if params["seq_mode"] else '*_{1,2}.fastq.gz'
    params["reads"] = os.path.join(relative_fastq_dir, glob_pattern)

    genome_data = parse_genome_file('genome.txt', current_dir)
    if not genome_data:
        sys.exit(1)
        
    selected_species = params["species"]
    selected_version = params["genome_version"]
    
    print("\n--- Preparing Genome Files ---")
    version_data = genome_data.get(selected_species, {}).get(selected_version, {})
    source = version_data.get('source', 'remote')

    if source == 'local':
        print(f"Using local genome for {selected_species} {selected_version}")
        params["fasta"] = version_data.get("fasta")
        params["gtf"] = version_data.get("gtf")
    else:
        fasta_url = version_data.get("fasta")
        gtf_url = version_data.get("gtf")

        genome_base_dir = os.path.join(current_dir, 'data', 'genomes', selected_species, selected_version)
        os.makedirs(genome_base_dir, exist_ok=True)

        if fasta_url:
            params["fasta"] = download_and_unzip(fasta_url, genome_base_dir)
        else:
            print(f"Warning: FASTA URL not found for {selected_species} {selected_version}", file=sys.stderr)
        
        if gtf_url:
            params["gtf"] = download_and_unzip(gtf_url, genome_base_dir)
        else:
            print(f"Warning: GTF URL not found for {selected_species} {selected_version}", file=sys.stderr)
    
    print("--- Genome Files Ready ---")

    print("\nFinal parameters for Nextflow:")
    print(json.dumps(params, indent=4))

    log_dir = os.path.join(current_dir, '.nextflow', 'log')
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
        "run",
        "main.nf",
        "-profile", "inside_container,standard",
        "-resume"
    ] + nextflow_params

    print("\nExecuting command:")
    print(" ".join(command))
    print("-" * 30)

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nError: Pipeline execution failed with exit code {e.returncode}")
    except FileNotFoundError:
        print("\nError: 'nextflow' command not found. Make sure you are in the correct environment.")

def main():
    """Main function to prepare parameters and launch the pipeline."""
    print("--- NGS Pipeline Launcher ---")
    args = get_args()
    current_dir = os.path.dirname(os.path.abspath(__file__))

    print("\nSetting up necessary directories...")
    if not setup_directories(current_dir):
        print("Failed to set up directories. Exiting.", file=sys.stderr)
        sys.exit(1)
    print("Directories setup successfully.")

    is_interactive = not (args.srr_list or args.fastq_dir or args.species or args.genome_version)

    if is_interactive:
        params, download_srr = handle_interactive_mode(current_dir)
    else:
        params, download_srr = handle_non_interactive_mode(args)

    confirm_and_run(params, download_srr, is_interactive, current_dir)

if __name__ == "__main__":
    main()
