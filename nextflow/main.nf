#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- 1. Parameter Definitions ---
params.reads = null
params.outdir = './results' // Set a default output directory
params.srr_list = null // This will be the primary input mechanism

// Genome-related parameters - MUST be provided by the user
params.fasta = null
params.gtf = null
params.species = null
params.genome_version = null

// Parameters for downstream analysis
params.run_de_analysis = false
params.meta_file = null
params.control_group = null
params.experiment_group = null


// --- 2. Validate Parameters & Create Input Channel ---

// Ensure essential genome parameters are provided
if (!params.fasta || !params.gtf || !params.species || !params.genome_version) {
    error """
    Missing essential genome parameters! Please provide all of the following:
    --fasta [path_to_fasta]
    --gtf [path_to_gtf]
    --species [species_name]
    --genome_version [version_name]
    """
}

// Create input channel from SRR list file
if (params.srr_list) {
    Channel
        .fromPath(params.srr_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { srr_id ->
            // Assuming a standard data directory structure
            def read1 = file("data/fastq/${srr_id}_1.fastq.gz", checkIfExists: true)
            def read2 = file("data/fastq/${srr_id}_2.fastq.gz", checkIfExists: true)
            tuple(srr_id, [read1, read2])
        }
        .ifEmpty { error "Cannot find any reads for SRR IDs in file: ${params.srr_list}" }
        .set { ch_reads }
} else if (params.reads) {
    // Fallback for manually providing reads pattern
    Channel
        .fromFilePairs(params.reads, size: 2, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { ch_reads }
} else {
    error "Input data not specified. Please provide --srr_list (a file with one SRR ID per line) or --reads (a glob pattern)."
}


// --- 3. Main Workflow Definition ---
workflow {
    // Define the expected path for the STAR index
    def star_index_path = file("${params.outdir}/genomes/${params.species}/${params.genome_version}/star_index")

    ch_star_index = Channel.empty()

    // Check if the index already exists
    if (star_index_path.exists() && star_index_path.isDirectory()) {
        log.info "Found existing STAR index at: ${star_index_path}"
        ch_star_index = Channel.fromPath(star_index_path)
    } else {
        log.info "STAR index not found. Generating a new one."
        // Step 1: Build STAR index
        STAR_GENOME_GENERATE(
            Channel.fromPath(params.fasta),
            Channel.fromPath(params.gtf)
        )
        ch_star_index = STAR_GENOME_GENERATE.out
    }
    
    // Step 2: Quality control and filtering with fastp
    FASTP(ch_reads)

    // Step 3: Combine fastp output with STAR index to ensure each sample gets the index path
    ch_for_alignment = FASTP.out.reads.combine(ch_star_index)
    STAR_ALIGN(ch_for_alignment)

    // Step 4: Collect all BAM files and perform unified quantification
    ch_bam_files = STAR_ALIGN.out.bam_ch.map { it[1] }.collect()
    FEATURECOUNTS(ch_bam_files, Channel.fromPath(params.gtf))

    // --- Optional Step 5: Downstream Differential Expression Analysis ---
    if (params.run_de_analysis) {
        if (!params.meta_file || !params.control_group || !params.experiment_group) {
            error "To run DE analysis, --meta_file, --control_group, and --experiment_group must be provided."
        }
        // Placeholder for future DE analysis module
    }
}

// --- 4. Process Definitions ---

process STAR_GENOME_GENERATE {
    publishDir "${params.outdir}/genomes/${params.species}/${params.genome_version}", mode: 'copy'
    
    label 'large_mem_process'

    input:
    path fasta
    path gtf

    output:
    path "star_index"

    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate \\
         --genomeDir star_index \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus}
    """
}

process FASTP {
    publishDir "${params.outdir}/fastp/${sample_id}", mode: 'copy'

    label 'default_process'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz"), emit: reads
    path "${sample_id}.html", emit: html
    path "${sample_id}.json", emit: json

    script:
    def (r1, r2) = reads
    def r1_out = "${sample_id}_trimmed_1.fastq.gz"
    def r2_out = r2 ? "${sample_id}_trimmed_2.fastq.gz" : ''
    def html_out = "${sample_id}.html"
    def json_out = "${sample_id}.json"
    def extra_args = r2 ? "-I ${r2} -O ${r2_out}" : ""

    """
    fastp -i ${r1} -o ${r1_out} -h ${html_out} -j ${json_out} --qualified_quality_phred 20 --unqualified_percent_limit 40 --length_required 36 ${extra_args}
    """
}

process STAR_ALIGN {
    publishDir "${params.outdir}/bam/${sample_id}", mode: 'copy', pattern: "*.{bam,log}"
    
    label 'large_mem_process'

    input:
    tuple val(sample_id), path(fastq_files), path(star_index)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam_ch
    path "*.log", emit: log_ch

    script:
    """
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${fastq_files.join(' ')} \\
         --readFilesCommand zcat \\
         --runThreadN ${task.cpus} \\
         --outFileNamePrefix "${sample_id}_" \\
         --outSAMtype BAM SortedByCoordinate
         
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    mv ${sample_id}_Log.final.out ${sample_id}.log
    """
}

process FEATURECOUNTS {
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    label 'default_process'

    input:
    path bams
    path gtf

    output:
    path "counts.txt"
    path "counts.txt.summary"

    script:
    def extra_params = "-p -s 2"
    """
    featureCounts -a ${gtf} \\
                  -o counts.txt \\
                  -T ${task.cpus} \\
                  ${extra_params} \\
                  ${bams.join(' ')}
    """
}
