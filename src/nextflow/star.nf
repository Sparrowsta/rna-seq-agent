#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * STAR RNA-seq比对流水线
 * 输入：FastP修剪后的FASTQ文件
 * 输出：比对后的BAM文件和统计信息
 */

// 参数定义 - 直接从 FastP 结果接收样本信息
params.sample_inputs = ""  // JSON格式的样本和FASTQ文件映射
params.star_index = ""
params.results_dir = ""
params.work_dir = ""

// STAR比对参数
params.outFilterMultimapNmax = 20
params.outFilterMismatchNoverReadLmax = 0.04
params.outFilterScoreMinOverLread = 0.33
params.outFilterMatchNminOverLread = 0.33
params.quantMode = "TranscriptomeSAM GeneCounts"
params.outSAMtype = "BAM SortedByCoordinate"
params.outSAMstrandField = "intronMotif"
params.twopassMode = "None"
params.limitBAMsortRAM = 4000000000

// 补充常用可调参数（与工具层对齐）
params.outSAMunmapped            = params.outSAMunmapped ?: "Within"
params.outSAMattributes          = params.outSAMattributes ?: "All"
params.outBAMsortingThreadN      = (params.outBAMsortingThreadN != null ? params.outBAMsortingThreadN as int : 1)
params.genomeLoad                = params.genomeLoad ?: "NoSharedMemory"
params.alignSJoverhangMin        = (params.alignSJoverhangMin != null ? params.alignSJoverhangMin as int : 8)
params.alignSJDBoverhangMin      = (params.alignSJDBoverhangMin != null ? params.alignSJDBoverhangMin as int : 1)
params.outFilterMismatchNmax     = (params.outFilterMismatchNmax != null ? params.outFilterMismatchNmax as int : 999)
params.alignIntronMin            = (params.alignIntronMin != null ? params.alignIntronMin as int : 20)
params.alignIntronMax            = (params.alignIntronMax != null ? params.alignIntronMax as int : 1000000)
params.alignMatesGapMax          = (params.alignMatesGapMax != null ? params.alignMatesGapMax as int : 1000000)

// 输入验证
if (!params.sample_inputs) error "Missing required parameter: sample_inputs"
if (!params.star_index) error "Missing required parameter: star_index"
if (!params.results_dir) error "Missing required parameter: results_dir"


// STAR比对进程
process STAR_ALIGN {
    cpus (params.resources?.star?.cpus ?: 8)
    memory (params.resources?.star?.memory ?: '32 GB')
    publishDir "${params.results_dir}/star", mode: 'copy'
    
    input:
    tuple val(sample_id), val(reads_info)
    path star_index
    
    output:
    tuple val(sample_id), path("${sample_id}/"), emit: aligned
    tuple val(sample_id), path("${sample_id}/${sample_id}.Log.final.out"), emit: logs
    
    script:
    def mem_gb = (params.limitBAMsortRAM / 1000000000).intValue()
    def is_paired = reads_info.is_paired
    def read1 = reads_info.read1
    def read2 = is_paired ? reads_info.read2 : ""
    
    """
    mkdir -p ${sample_id}

    # 构建STAR命令，根据单端/双端数据调整（不进入子目录）
    micromamba run -n align_env STAR \
        --genomeDir ${star_index} \
        --readFilesIn ${read1}${is_paired ? ' ' + read2 : ''} \
        --readFilesCommand ${params.readFilesCommand ?: 'zcat'} \
        --runThreadN ${task.cpus} \
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \
        --outFilterMismatchNoverReadLmax ${params.outFilterMismatchNoverReadLmax} \
        --outFilterMismatchNmax ${params.outFilterMismatchNmax} \
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \
        --alignSJoverhangMin ${params.alignSJoverhangMin} \
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \
        --alignIntronMin ${params.alignIntronMin} \
        --alignIntronMax ${params.alignIntronMax} \
        --alignMatesGapMax ${params.alignMatesGapMax} \
        --quantMode ${params.quantMode} \
        --outSAMtype ${params.outSAMtype} \
        --outSAMstrandField ${params.outSAMstrandField} \
        --outSAMunmapped ${params.outSAMunmapped} \
        --outSAMattributes ${params.outSAMattributes} \
        --twopassMode ${params.twopassMode} \
        --limitBAMsortRAM ${params.limitBAMsortRAM} \
        --outBAMsortingThreadN ${params.outBAMsortingThreadN} \
        --genomeLoad ${params.genomeLoad} \
        --outFileNamePrefix ${sample_id}/${sample_id}.
    """
}

// 主工作流
workflow {
    // 直接从参数解析样本输入信息
    def sample_data = new groovy.json.JsonSlurper().parseText(params.sample_inputs)
    
    // 创建样本-读取信息的输入通道
    input_files = Channel.from(sample_data)
        .map { sample_info ->
            def sample_id = sample_info.sample_id
            def reads_info = [
                is_paired: sample_info.is_paired,
                read1: file(sample_info.read1),
                read2: sample_info.is_paired ? file(sample_info.read2) : null
            ]
            return tuple(sample_id, reads_info)
        }
    
    // 执行STAR比对（将索引作为常量广播到所有样本）
    star_index_ch = Channel.value(file(params.star_index))
    STAR_ALIGN(input_files, star_index_ch)
    
    // 收集日志信息
    STAR_ALIGN.out.logs.collectFile(name: 'star_summary.txt', newLine: true)
}
