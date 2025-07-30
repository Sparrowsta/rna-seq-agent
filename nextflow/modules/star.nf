//
// STAR Alignment Workflow
// Responsibility: Perform sequence alignment using STAR
//
nextflow.enable.dsl=2

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

workflow STAR_WORKFLOW {
    take:
        reads_ch
        star_index_ch
    
    main:
        STAR_ALIGN(reads_ch.combine(star_index_ch))
    
    emit:
        bam = STAR_ALIGN.out.bam_ch
        log = STAR_ALIGN.out.log_ch
}

// 如果直接运行此脚本
if (params.reads && params.star_index) {
    Channel
        .fromFilePairs(params.reads, size: -1)
        .map { sample_id, files -> [sample_id, files.sort()] }
        .set { reads_ch }
    
    Channel
        .fromPath(params.star_index)
        .set { star_index_ch }
    
    STAR_WORKFLOW(reads_ch, star_index_ch)
} 