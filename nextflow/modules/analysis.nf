//
// Module: Core RNA-Seq Analysis
// Responsibility: Perform QC, Alignment, and Quantification.
// Assumes it receives only the samples that need to be processed.
//
nextflow.enable.dsl=2

workflow ANALYSIS {
    take:
        reads_ch      // channel: [ [id, [/path/to/r1, /path/to/r2]], ... ]
        star_index_ch // channel: /path/to/star_index
        gtf_ch        // channel: /path/to/genes.gtf

    main:
        // 1. Quality Control
        FASTP(reads_ch)

        // 2. Alignment
        ch_for_alignment = FASTP.out.reads.combine(star_index_ch)
        STAR_ALIGN(ch_for_alignment)

        // 3. Quantification
        ch_bam_files = STAR_ALIGN.out.bam_ch.map { it[1] }.collect()
        FEATURECOUNTS(ch_bam_files, gtf_ch)

    emit:
        counts = FEATURECOUNTS.out.counts
}


// --- Process Definitions ---

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
    path "counts.txt", emit: counts
    path "counts.txt.summary", emit: summary

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