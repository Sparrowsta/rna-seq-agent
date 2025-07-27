#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- 1. Parameter Definitions ---
params.reads = null
params.outdir = './results'      // Set a default output directory

// Genome-related parameters - MUST be provided by the agent
params.fasta = null
params.gtf = null
params.species = null
params.genome_version = null

// Parameters for tools
params.fastp_q_val = 20
params.fastp_unqual_pct = 40
params.fastp_len_req = 50
params.fc_is_paired_end = true
params.fc_strand_spec = 2

// Parameters for downstream analysis
params.run_de_analysis = false
params.meta_file = null
params.control_group = null
params.experiment_group = null


// --- 2. Workflow Definition ---
workflow {
    // --- 1. Input Channel Setup ---
    if (!params.reads) {
        error "Input reads not specified. Please provide a glob pattern for FASTQ files via --reads, e.g., 'path/to/*_{1,2}.fastq.gz'"
    }

    ch_reads = Channel
        .fromFilePairs(params.reads, size: -1) { file ->
            file.name.replaceAll(/_?[12]\.fastq\.gz$|_\.fastq\.gz$/, "")
        }
        .map { sample_id, files ->
            tuple(sample_id, files.sort())
        }

    // --- 2. Genome Index Management ---
    def published_index = file(params.fasta).parent.resolve('star_index')
    ch_star_index = Channel.empty()
    if (published_index.exists() && published_index.isDirectory()) {
        log.info "Found published STAR index: ${published_index}"
        ch_star_index = Channel.fromPath(published_index)
    } else {
        log.info "Published STAR index not found, generating..."
        STAR_GENOME_GENERATE(
            Channel.fromPath(params.fasta),
            Channel.fromPath(params.gtf)
        )
        ch_star_index = STAR_GENOME_GENERATE.out
    }

    // --- 3. 智能 Quality Control ---
    ch_reads.filter { sid, _ ->
        !file("${params.outdir}/fastp/${sid}/${sid}_trimmed_1.fastq.gz").exists()
    }.set { ch_reads_to_run_fastp }

    ch_reads.filter { sid, _ ->
        file("${params.outdir}/fastp/${sid}/${sid}_trimmed_1.fastq.gz").exists()
    }.set { ch_reads_existing_fastp }

    FASTP(ch_reads_to_run_fastp)

    ch_fastp_results_from_existing = ch_reads_existing_fastp.map { sid, _ ->
        log.info "Skipping FASTP for ${sid}: published files found."
        def r1 = file("${params.outdir}/fastp/${sid}/${sid}_trimmed_1.fastq.gz")
        def r2 = file("${params.outdir}/fastp/${sid}/${sid}_trimmed_2.fastq.gz")
        return tuple(sid, r2.exists() ? [r1, r2] : [r1])
    }
    ch_trimmed_reads = FASTP.out.reads.mix(ch_fastp_results_from_existing)

    // --- 4. 智能 Alignment ---
    ch_for_alignment = ch_trimmed_reads.combine(ch_star_index)

    ch_for_alignment.filter { sid, _reads, _index ->
        !file("${params.outdir}/bam/${sid}/${sid}.bam").exists()
    }.set { ch_to_run_align }

    ch_for_alignment.filter { sid, _reads, _index ->
        file("${params.outdir}/bam/${sid}/${sid}.bam").exists()
    }.set { ch_existing_align }

    STAR_ALIGN(ch_to_run_align)

    ch_star_results_from_existing = ch_existing_align.map { sid, _reads, _index ->
        log.info "Skipping STAR Alignment for ${sid}: published BAM found."
        return tuple(sid, file("${params.outdir}/bam/${sid}/${sid}.bam"))
    }
    ch_aligned_bams = STAR_ALIGN.out.bam_ch.mix(ch_star_results_from_existing)

    // --- 5. 智能 Quantification ---
    def final_counts_file = file("${params.outdir}/featurecounts/counts.txt")
    if (final_counts_file.exists()) {
        log.info "Final counts file found. Skipping quantification."
    } else {
        // .collect() is necessary here as featureCounts takes a list of files
        ch_bam_files = ch_aligned_bams.map { it[1] }.collect()
        FEATURECOUNTS(ch_bam_files, Channel.fromPath(params.gtf))
    }

    // --- 6. (Optional) Differential Expression Analysis ---
    if (params.run_de_analysis) {
        if (!params.meta_file || !params.control_group || !params.experiment_group) {
            error "To run DE analysis, --meta_file, --control_group, and --experiment_group must be provided."
        }
        // Placeholder for future DE analysis module
    }
}

// --- 3. Process Definitions ---

process STAR_GENOME_GENERATE {
    publishDir path: "${file(params.fasta).parent}", mode: 'copy'
    
    label 'large_mem_process'

    input:
    path fasta
    path gtf

    output:
    path "star_index"

    script:
    """
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
    source activate qc_env
    fastp -i ${r1} -o ${r1_out} -h ${html_out} -j ${json_out} --qualified_quality_phred ${params.fastp_q_val} --unqualified_percent_limit ${params.fastp_unqual_pct} --length_required ${params.fastp_len_req} ${extra_args}
    """
}

process STAR_ALIGN {
    publishDir "${params.outdir}/bam/${sample_id}", mode: 'copy', pattern: "*.{bam,log}"
    
    label 'large_mem_process'

    input:
    tuple val(sample_id), path(reads), path(star_index)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_ch
    path "${sample_id}.log", emit: log_ch

    script:
    def read_files_str = reads instanceof List ? reads.join(' ') : reads
    """
    source activate align_env
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${read_files_str} \\
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
    def extra_params = ""
    if (params.fc_is_paired_end) {
        extra_params += " -p"
    }
    extra_params += " -s ${params.fc_strand_spec}"
    """
    source activate quant_env
    featureCounts -a ${gtf} \\
                  -o counts.txt \\
                  -T ${task.cpus} \\
                  ${extra_params} \\
                  ${bams.join(' ')}
    """
}
