//
// Module: Download Genome Source Files
// Responsibility: Unconditionally download FASTA and GTF files from URLs.
//
nextflow.enable.dsl=2

workflow DOWNLOAD_GENOME {
    take:
        genome_info_ch // channel: [ name, {species:..., fasta_url:..., gtf_url:...} ]

    main:
        WGET_GENOME_FILES(genome_info_ch)

    emit:
        // It emits the paths to the downloaded FASTA and GTF files,
        // ensuring downstream processes can depend on them directly.
        genome_files = WGET_GENOME_FILES.out
}

process WGET_GENOME_FILES {
    // Publish downloaded files to the genomes directory
    publishDir "data/genomes/${params.species}/${params.genome_version}/", mode: 'copy'
    
    label 'default_process'

    input:
    tuple val(genome_name), val(genome_info)

    output:
    tuple val(genome_name), path(fasta_path), path(gtf_path)

    script:
    // Define the output paths based on the config. These variables are now
    // available to the 'output' block.
    // Calculate the final paths (after decompression)
    // Use simple filenames to avoid path duplication
    fasta_path = genome_info.fasta.split('/')[-1]  // Just the filename
    gtf_path = genome_info.gtf.split('/')[-1]      // Just the filename
    // Define temporary paths for compressed files
    temp_fasta = fasta_path + ".gz"
    temp_gtf = gtf_path + ".gz"
    """
    #!/bin/bash
    set -e
    
    echo "Preparing directories for ${genome_name}..."
    mkdir -p \$(dirname ${fasta_path})
    mkdir -p \$(dirname ${gtf_path})

    echo "Downloading FASTA from ${genome_info.fasta_url}..."
    # Download to temporary .gz file first
    wget --quiet -O ${temp_fasta} "${genome_info.fasta_url}"

    echo "Downloading GTF from ${genome_info.gtf_url}..."
    # Download to temporary .gz file first
    wget --quiet -O ${temp_gtf} "${genome_info.gtf_url}"

    echo "Decompressing files..."
    # Decompress FASTA
    echo "Decompressing FASTA file..."
    gunzip -f "${temp_fasta}"
    
    # Decompress GTF
    echo "Decompressing GTF file..."
    gunzip -f "${temp_gtf}"

    echo "Downloads for ${genome_name} complete."
    """
}