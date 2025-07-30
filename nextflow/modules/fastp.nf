//
// FASTP Quality Control Workflow
// Responsibility: Perform quality control on FASTQ files
//
nextflow.enable.dsl=2

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

workflow FASTP_WORKFLOW {
    take:
        reads_ch
    
    main:
        FASTP(reads_ch)
    
    emit:
        reads = FASTP.out.reads
        html = FASTP.out.html
        json = FASTP.out.json
}

// 如果直接运行此脚本
if (params.reads) {
    Channel
        .fromFilePairs(params.reads, size: -1)
        .map { sample_id, files -> [sample_id, files.sort()] }
        .set { reads_ch }
    
    FASTP_WORKFLOW(reads_ch)
} 