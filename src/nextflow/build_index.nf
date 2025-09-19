#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * STAR索引构建流水线
 * 输入：基因组FASTA和GTF文件
 * 输出：STAR索引目录
 */

// 参数定义
params.genome_fasta = ""
params.genome_gtf = ""
params.star_index_dir = ""
params.sjdb_overhang = 100
params.limitGenomeGenerateRAM = 32000000000

// 输入验证
if (!params.genome_fasta) error "Missing required parameter: genome_fasta"
if (!params.genome_gtf) error "Missing required parameter: genome_gtf"
if (!params.star_index_dir) error "Missing required parameter: star_index_dir"


// STAR索引构建进程
process BUILD_STAR_INDEX {
    cpus (params.resources?.star?.cpus ?: 8)
    memory (params.resources?.star?.memory ?: '32 GB')
    publishDir "${params.star_index_dir}", mode: 'copy'
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "*", emit: index_files
    
    script:
    """
    # 获取输入文件的完整路径，避免链接冲突
    FASTA_FILE=\$(readlink -f ${genome_fasta})
    GTF_FILE=\$(readlink -f ${genome_gtf})

    
    # STAR索引构建 - 直接在当前目录生成索引文件  
    micromamba run -n align_env STAR \\
        --runMode genomeGenerate \\
        --genomeDir . \\
        --genomeFastaFiles \$FASTA_FILE \\
        --sjdbGTFfile \$GTF_FILE \\
        --sjdbOverhang ${params.sjdb_overhang} \\
        --runThreadN ${task.cpus} \\
        --limitGenomeGenerateRAM ${params.limitGenomeGenerateRAM}
    
    # 创建完成标记
    echo "STAR index built at \$(date)" > Log.out
    """
}

// 主工作流
workflow {
    // 输入文件通道
    genome_fasta_ch = Channel.fromPath(params.genome_fasta)
    genome_gtf_ch = Channel.fromPath(params.genome_gtf)
    
    // 构建索引
    BUILD_STAR_INDEX(genome_fasta_ch, genome_gtf_ch)
}
