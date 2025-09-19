#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * HISAT2 RNA-seq比对流水线
 * 输入：FastP修剪后的FASTQ文件
 * 输出：比对后的BAM文件和统计信息
 */

// 参数定义 - 直接从 FastP 结果接收样本信息
params.sample_inputs = ""  // JSON格式的样本和FASTQ文件映射
params.hisat2_index = ""   // HISAT2索引目录前缀
params.results_dir = ""
params.work_dir = ""

// HISAT2比对基础参数
params.rna_strandness = "unstranded"        // 链特异性: unstranded/RF/FR
params.dta = true                           // 启用下游转录本组装模式
// 线程数由 process 指令 cpus 管理，命令行使用 ${task.cpus}
params.max_intronlen = 500000               // 最大内含子长度
params.min_intronlen = 20                   // 最小内含子长度

// HISAT2比对质量控制参数（与 DEFAULT_HISAT2_PARAMS 对齐）
params.score_min = params.score_min ?: "L,0,-0.2"      // 最小比对得分函数
params.mp = params.mp ?: "6,2"                          // 错配惩罚: max,min
params.sp = params.sp ?: "2,1"                          // 软剪切惩罚: max,min
params.np = (params.np != null ? params.np as int : 1)  // N碱基惩罚
params.rdg = params.rdg ?: "5,3"                        // read gap开放和延伸惩罚
params.rfg = params.rfg ?: "5,3"                        // 参考基因组gap开放和延伸惩罚

// 比对控制参数
params.no_mixed = (params.no_mixed != null ? params.no_mixed as boolean : false)        // 允许混合比对
params.no_discordant = (params.no_discordant != null ? params.no_discordant as boolean : false)  // 允许不一致比对
params.k = (params.k != null ? params.k as int : 5)                                     // 最多报告k个比对结果
params.n_ceil = params.n_ceil ?: "L,0,0.15"             // N碱基数量上限函数

// 种子和扩展参数
params.i = params.i ?: "S,1,1.15"                       // 种子间隔函数
params.L = (params.L != null ? params.L as int : 20)    // 种子长度
params.N = (params.N != null ? params.N as int : 0)     // 每个种子最多错配数
params.D = (params.D != null ? params.D as int : 15)    // 扩展尝试次数
params.R = (params.R != null ? params.R as int : 2)     // 重新播种次数

// 输出控制参数
params.no_unal = (params.no_unal != null ? params.no_unal as boolean : false)  // 不输出未比对reads
params.quiet = (params.quiet != null ? params.quiet as boolean : false)        // 安静模式

// 输入验证
if (!params.sample_inputs) error "Missing required parameter: sample_inputs"
if (!params.hisat2_index) error "Missing required parameter: hisat2_index"
if (!params.results_dir) error "Missing required parameter: results_dir"

// HISAT2比对进程
process HISAT2_ALIGN {
    cpus (params.resources?.hisat2?.cpus ?: 8)
    memory (params.resources?.hisat2?.memory ?: '32 GB')
    publishDir "${params.results_dir}/hisat2", mode: 'copy'
    
    input:
    tuple val(sample_id), val(reads_info)
    val hisat2_index_prefix
    
    output:
    tuple val(sample_id), path("${sample_id}/"), emit: aligned
    tuple val(sample_id), path("${sample_id}/${sample_id}.align_summary.txt"), emit: logs
    
    script:
    def is_paired = reads_info.is_paired
    def read1 = reads_info.read1
    def read2 = is_paired ? reads_info.read2 : ""
    
    // 构建HISAT2基础命令参数
    def hisat2_cmd = "micromamba run -n align_env hisat2"
    def basic_params = [
        "-x ${hisat2_index_prefix}",
        "-p ${task.cpus}",
        "--max-intronlen ${params.max_intronlen}",
        "--min-intronlen ${params.min_intronlen}"
    ]
    
    // RNA-seq 特定参数
    if (params.rna_strandness && params.rna_strandness != "unstranded") {
        basic_params.add("--rna-strandness ${params.rna_strandness}")
    }
    if (params.dta) {
        basic_params.add("--dta")
    }
    
    // 质量控制参数
    def quality_params = []
    if (params.score_min) quality_params.add("--score-min ${params.score_min}")
    if (params.mp) quality_params.add("--mp ${params.mp}")
    if (params.sp) quality_params.add("--sp ${params.sp}")
    if (params.np != null) quality_params.add("--np ${params.np}")
    if (params.rdg) quality_params.add("--rdg ${params.rdg}")
    if (params.rfg) quality_params.add("--rfg ${params.rfg}")
    
    // 比对控制参数
    def control_params = []
    if (params.no_mixed) control_params.add("--no-mixed")
    if (params.no_discordant) control_params.add("--no-discordant")
    if (params.k) control_params.add("-k ${params.k}")
    if (params.n_ceil) control_params.add("--n-ceil ${params.n_ceil}")
    
    // 种子参数
    def seed_params = []
    if (params.i) seed_params.add("-i ${params.i}")
    if (params.L) seed_params.add("-L ${params.L}")
    if (params.N) seed_params.add("-N ${params.N}")
    if (params.D) seed_params.add("-D ${params.D}")
    if (params.R) seed_params.add("-R ${params.R}")
    
    // 输出控制参数
    def output_params = []
    if (params.no_unal) output_params.add("--no-unal")
    if (params.quiet) output_params.add("--quiet")
    
    // 组合所有参数
    def all_params = basic_params + quality_params + control_params + seed_params + output_params
    
    """
    mkdir -p ${sample_id}
    
    # 构建HISAT2命令，根据单端/双端数据调整
    ${hisat2_cmd} \\
        ${all_params.join(' \\\n        ')} \\
        ${is_paired ? "-1 ${read1} -2 ${read2}" : "-U ${read1}"} \\
        -S ${sample_id}/${sample_id}.hisat2.sam \\
        2> ${sample_id}/${sample_id}.align_summary.txt
    
    # 转换SAM为BAM并排序
    micromamba run -n align_env samtools view -@ ${task.cpus} -bS ${sample_id}/${sample_id}.hisat2.sam | \\
    micromamba run -n align_env samtools sort -@ ${task.cpus} -o ${sample_id}/${sample_id}.hisat2.bam
    
    # 删除中间SAM文件节省空间
    rm ${sample_id}/${sample_id}.hisat2.sam
    
    # 生成BAM索引
    micromamba run -n align_env samtools index ${sample_id}/${sample_id}.hisat2.bam
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
    
    // 执行HISAT2比对（将索引前缀作为常量广播到所有样本）
    hisat2_index_ch = Channel.value(params.hisat2_index)
    HISAT2_ALIGN(input_files, hisat2_index_ch)
    
    // 收集日志信息
    HISAT2_ALIGN.out.logs.collectFile(name: 'hisat2_summary.txt', newLine: true)
}
