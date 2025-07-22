#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- 1. Parameter Definitions ---
params.reads = null
params.outdir = 'data'
params.seq_mode = false

// Genome-related parameters
params.fasta = null
params.gtf = null
params.species = null
params.genome_version = null

// --- 2. Create Input Channel ---
if (params.seq_mode) {
    Channel
        .fromFilePairs(params.reads, size: 1, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { ch_reads }
} else {
    Channel
        .fromFilePairs(params.reads, size: 2, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { ch_reads }
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
    FEATURECOUNTS(ch_bam_files)
}

// --- 4. Process Definitions ---

process STAR_GENOME_GENERATE {
    // The index is large, so we publish it to the genome's directory
    // to keep related files together.
    publishDir "${params.outdir}/genomes/${params.species}/${params.genome_version}", mode: 'copy'
    
    label 'large_mem_process'

    input:
    path fasta
    path gtf

    output:
    path "star_index"

    script:
    """
    echo "Activating align_env to build STAR index..."
    source activate align_env

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
    source activate qc_env
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
    echo "Activating align_env for STAR on sample ${sample_id}..."
    source activate align_env
    
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${fastq_files.join(' ')} \\
         --readFilesCommand zcat \\
         --runThreadN ${task.cpus} \\
         --outFileNamePrefix "${sample_id}_" \\
         --outSAMtype BAM SortedByCoordinate
         
    # Rename output files to remove the extra suffix for consistency
    mv ${sample_id}_Aligned.sortedByCoord.out.bam ${sample_id}.bam
    mv ${sample_id}_Log.final.out ${sample_id}.log
    """
}

process FEATURECOUNTS {
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    val bams

    output:
    path "counts.txt"
    path "counts.txt.summary"

    script:
    def extra_params = params.seq_mode ? "-s 0" : "-p -s 2"
    """
    echo "Activating quant_env for featureCounts..."
    source activate quant_env
    
    featureCounts -a ${params.gtf} \\
                  -o counts.txt \\
                  -T ${task.cpus} \\
                  ${extra_params} \\
                  ${bams.join(' ')}
    """
}