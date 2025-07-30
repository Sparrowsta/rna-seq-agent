//
// FeatureCounts Quantification Workflow
// Responsibility: Perform gene quantification using featureCounts
//
nextflow.enable.dsl=2

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
    
    // 确保bams是一个列表
    def bam_files = bams instanceof List ? bams : [bams]
    def bam_args = bam_files.join(' ')
    
    """
    source activate quant_env
    featureCounts -a ${gtf} \\
                  -o counts.txt \\
                  -T ${task.cpus} \\
                  ${extra_params} \\
                  ${bam_args}
    """
}

workflow FEATURECOUNTS_WORKFLOW {
    take:
        bams_ch
        gtf_ch
    
    main:
        FEATURECOUNTS(bams_ch, gtf_ch)
    
    emit:
        counts = FEATURECOUNTS.out.counts
        summary = FEATURECOUNTS.out.summary
}

// 如果直接运行此脚本
if (params.bams && params.gtf) {
    Channel
        .fromPath(params.bams.split(','))
        .set { bams_ch }
    
    Channel
        .fromPath(params.gtf)
        .set { gtf_ch }
    
    FEATURECOUNTS_WORKFLOW(bams_ch, gtf_ch)
} 