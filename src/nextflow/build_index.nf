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
params.runThreadN = 4
params.limitGenomeGenerateRAM = 32000000000

// 输入验证
if (!params.genome_fasta) error "Missing required parameter: genome_fasta"
if (!params.genome_gtf) error "Missing required parameter: genome_gtf"
if (!params.star_index_dir) error "Missing required parameter: star_index_dir"

log.info """
=================================
STAR 索引构建流水线
=================================
基因组FASTA: ${params.genome_fasta}
基因组GTF: ${params.genome_gtf}
索引输出目录: ${params.star_index_dir}
sjdbOverhang: ${params.sjdb_overhang}
=================================
"""

// STAR索引构建进程
process BUILD_STAR_INDEX {
    publishDir "${params.star_index_dir}", mode: 'copy'
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "*", emit: index_files
    
    script:
    """
    mkdir -p star_index_temp
    
    micromamba run -n align_env STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${params.runThreadN} \\
        --genomeDir star_index_temp \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${genome_gtf} \\
        --sjdbOverhang ${params.sjdb_overhang} \\
        --limitGenomeGenerateRAM ${params.limitGenomeGenerateRAM}
    
    # 移动索引文件到当前目录
    mv star_index_temp/* .
    rmdir star_index_temp
    
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
