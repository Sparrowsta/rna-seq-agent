#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * HISAT2索引构建流水线
 * 输入：基因组FASTA和GTF文件
 * 输出：HISAT2索引目录
 */

// 参数定义
params.genome_fasta = ""
params.genome_gtf = ""
params.hisat2_index_dir = ""
params.index_basename = "genome"       // 索引文件前缀名称
params.p = 4                           // 线程数
params.snp = ""                        // SNP文件路径（可选）
params.haplotype = ""                  // 单倍型文件路径（可选）
params.ss = ""                         // 剪接位点文件路径（可选，可从GTF生成）
params.exon = ""                       // 外显子文件路径（可选，可从GTF生成）

// 输入验证
if (!params.genome_fasta) error "Missing required parameter: genome_fasta"
if (!params.hisat2_index_dir) error "Missing required parameter: hisat2_index_dir"

// 从GTF提取剪接位点和外显子信息进程
process EXTRACT_SPLICE_SITES {
    input:
    path genome_gtf
    
    output:
    path "splice_sites.txt", emit: splice_sites, optional: true
    path "exons.txt", emit: exons, optional: true
    
    when:
    genome_gtf.name != 'NO_FILE'
    
    script:
    """
    # 检查GTF文件是否存在且不为空
    if [[ -s ${genome_gtf} ]]; then
        # 提取剪接位点
        micromamba run -n align_env hisat2_extract_splice_sites.py ${genome_gtf} > splice_sites.txt
        
        # 提取外显子
        micromamba run -n align_env hisat2_extract_exons.py ${genome_gtf} > exons.txt
        
        echo "Extracted splice sites and exons from GTF file"
    else
        echo "No GTF file provided or file is empty, skipping splice site extraction"
        touch splice_sites.txt
        touch exons.txt
    fi
    """
}

// HISAT2索引构建进程
process BUILD_HISAT2_INDEX {
    publishDir "${params.hisat2_index_dir}", mode: 'copy'
    
    input:
    path genome_fasta
    path splice_sites
    path exons
    
    output:
    path "*.ht2", emit: index_files
    path "build_hisat2_index.log", emit: build_log
    
    script:
    def splice_sites_param = splice_sites.name != 'NO_FILE' && splice_sites.size() > 0 ? "--ss ${splice_sites}" : ""
    def exons_param = exons.name != 'NO_FILE' && exons.size() > 0 ? "--exon ${exons}" : ""
    def snp_param = params.snp ? "--snp ${params.snp}" : ""
    def haplotype_param = params.haplotype ? "--haplotype ${params.haplotype}" : ""
    
    """
    # 获取输入文件的完整路径，避免链接冲突
    FASTA_FILE=\$(readlink -f ${genome_fasta})
    
    # 记录开始时间
    echo "Starting HISAT2 index build at \$(date)" > build_hisat2_index.log
    echo "Genome FASTA: \$FASTA_FILE" >> build_hisat2_index.log
    echo "Index basename: ${params.index_basename}" >> build_hisat2_index.log
    echo "Threads: ${params.p}" >> build_hisat2_index.log
    
    # 清理当前目录的输入文件软链接
    rm -f ${genome_fasta}
    if [[ -f ${splice_sites} && "${splice_sites.name}" != "NO_FILE" ]]; then
        echo "Using splice sites file: ${splice_sites}" >> build_hisat2_index.log
    fi
    if [[ -f ${exons} && "${exons.name}" != "NO_FILE" ]]; then
        echo "Using exons file: ${exons}" >> build_hisat2_index.log
    fi
    
    # HISAT2索引构建
    micromamba run -n align_env hisat2-build \\
        ${splice_sites_param} \\
        ${exons_param} \\
        ${snp_param} \\
        ${haplotype_param} \\
        -p ${params.p} \\
        \$FASTA_FILE \\
        ${params.index_basename} \\
        2>&1 | tee -a build_hisat2_index.log
    
    # 检查索引构建是否成功
    if ls ${params.index_basename}.*.ht2 1> /dev/null 2>&1; then
        echo "HISAT2 index built successfully at \$(date)" >> build_hisat2_index.log
        echo "Index files:" >> build_hisat2_index.log
        ls -la ${params.index_basename}.*.ht2 >> build_hisat2_index.log
    else
        echo "ERROR: HISAT2 index build failed!" >> build_hisat2_index.log
        exit 1
    fi
    """
}

// 主工作流
workflow {
    // 输入文件通道
    genome_fasta_ch = Channel.fromPath(params.genome_fasta)
    
    // GTF文件通道 - 如果未提供则创建空文件标记
    if (params.genome_gtf && params.genome_gtf != "") {
        genome_gtf_ch = Channel.fromPath(params.genome_gtf)
    } else {
        genome_gtf_ch = Channel.value(file('NO_FILE'))
    }
    
    // 提取剪接位点和外显子信息（如果有GTF文件）
    EXTRACT_SPLICE_SITES(genome_gtf_ch)
    
    // 构建索引 - 如果没有GTF，splice_sites和exons为空文件
    splice_sites_ch = EXTRACT_SPLICE_SITES.out.splice_sites.ifEmpty(file('NO_FILE'))
    exons_ch = EXTRACT_SPLICE_SITES.out.exons.ifEmpty(file('NO_FILE'))
    
    BUILD_HISAT2_INDEX(genome_fasta_ch, splice_sites_ch, exons_ch)
}