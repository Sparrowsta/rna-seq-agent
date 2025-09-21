nextflow.enable.dsl=2

/*
 * HISAT2索引构建流水线
 * 输入：基因组FASTA和GTF文件
 * 输出：HISAT2索引目录
 * 支持所有HISAT2索引构建参数的动态传递
 */

// 基础必需参数
params.genome_fasta = ""
params.genome_gtf = ""
params.hisat2_index_dir = ""
params.index_basename = "genome"       // 索引文件前缀名称

// HISAT2 索引构建参数 - 支持动态传递
params.p = 4                           // 线程数
params.runThreadN = 8                  // 兼容命名
params.large_index = false             // 大索引模式
params.offrate = 4                     // 偏移速率

// 剪接位点和外显子文件（可选）
params.ss = null                       // 剪接位点文件路径
params.exon = null                     // 外显子文件路径
params.snp = null                      // SNP文件路径（可选）
params.haplotype = null                // 单倍型文件路径（可选）

// 高级索引参数（可选）
params.ftabchars = null                // ftab字符数
params.local = false                   // 本地构建模式
params.packed = false                  // 压缩索引
params.bmax = null                     // 最大bucket大小
params.bmaxdivn = null                 // bmax/参考长度比值
params.dcv = null                      // 差异覆盖采样
params.nodc = false                    // 不使用差异覆盖
params.noref = false                   // 不构建参考字符串
params.justref = false                 // 仅构建参考字符串
params.seed = null                     // 随机种子
params.cutoff = null                   // 截止值

// 输入验证
if (!params.genome_fasta) error "Missing required parameter: genome_fasta"
if (!params.hisat2_index_dir) error "Missing required parameter: hisat2_index_dir"

// 从GTF提取剪接位点和外显子信息进程
process EXTRACT_SPLICE_SITES {
    cpus (params.resources?.hisat2?.cpus ?: 4)
    memory (params.resources?.hisat2?.memory ?: '8 GB')
    
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
        echo "GTF file is empty or does not exist, skipping splice site extraction"
        touch splice_sites.txt
        touch exons.txt
    fi
    """
}

// HISAT2索引构建进程
process BUILD_HISAT2_INDEX {
    cpus (params.resources?.hisat2?.cpus ?: params.p ?: params.runThreadN ?: 4)
    memory (params.resources?.hisat2?.memory ?: '16 GB')
    publishDir "${params.hisat2_index_dir}", mode: 'copy'
    
    input:
    path genome_fasta
    path splice_sites
    path exons
    
    output:
    path "*.ht2*", emit: index_files
    
    script:
    """
    # 获取输入文件的完整路径
    FASTA_FILE=\$(readlink -f ${genome_fasta})
    
    # 构建HISAT2索引命令 - 支持动态参数
    HISAT2_CMD="micromamba run -n align_env hisat2-build"
    
    # 添加基础参数
    HISAT2_CMD="\$HISAT2_CMD -p ${task.cpus}"
    HISAT2_CMD="\$HISAT2_CMD --offrate ${params.offrate}"
    
    # 添加大索引模式
    if [[ "${params.large_index}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --large-index"
    fi
    
    # 添加剪接位点文件（如果可用）
    if [[ -s splice_sites.txt ]]; then
        HISAT2_CMD="\$HISAT2_CMD --ss splice_sites.txt"
    elif [[ "${params.ss}" != "null" && "${params.ss}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --ss ${params.ss}"
    fi
    
    # 添加外显子文件（如果可用）
    if [[ -s exons.txt ]]; then
        HISAT2_CMD="\$HISAT2_CMD --exon exons.txt"
    elif [[ "${params.exon}" != "null" && "${params.exon}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --exon ${params.exon}"
    fi
    
    # 添加可选高级参数（仅在非null时添加）
    if [[ "${params.snp}" != "null" && "${params.snp}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --snp ${params.snp}"
    fi
    
    if [[ "${params.haplotype}" != "null" && "${params.haplotype}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --haplotype ${params.haplotype}"
    fi
    
    if [[ "${params.ftabchars}" != "null" && "${params.ftabchars}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --ftabchars ${params.ftabchars}"
    fi
    
    if [[ "${params.local}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --local"
    fi
    
    if [[ "${params.packed}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --packed"
    fi
    
    if [[ "${params.bmax}" != "null" && "${params.bmax}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --bmax ${params.bmax}"
    fi
    
    if [[ "${params.bmaxdivn}" != "null" && "${params.bmaxdivn}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --bmaxdivn ${params.bmaxdivn}"
    fi
    
    if [[ "${params.dcv}" != "null" && "${params.dcv}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --dcv ${params.dcv}"
    fi
    
    if [[ "${params.nodc}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --nodc"
    fi
    
    if [[ "${params.noref}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --noref"
    fi
    
    if [[ "${params.justref}" == "true" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --justref"
    fi
    
    if [[ "${params.seed}" != "null" && "${params.seed}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --seed ${params.seed}"
    fi
    
    if [[ "${params.cutoff}" != "null" && "${params.cutoff}" != "" ]]; then
        HISAT2_CMD="\$HISAT2_CMD --cutoff ${params.cutoff}"
    fi
    
    # 添加输入文件和输出前缀
    HISAT2_CMD="\$HISAT2_CMD \$FASTA_FILE ${params.index_basename}"
    
    # 执行HISAT2索引构建
    echo "执行HISAT2索引构建命令: \$HISAT2_CMD"
    eval \$HISAT2_CMD
    
    # 创建完成标记
    echo "HISAT2 index built at \$(date)" > ${params.index_basename}.log
    """
}

// 主工作流
workflow {
    // 输入文件通道
    genome_fasta_ch = Channel.fromPath(params.genome_fasta)
    
    // 处理GTF文件（可选）
    if (params.genome_gtf && params.genome_gtf != "") {
        genome_gtf_ch = Channel.fromPath(params.genome_gtf)
        EXTRACT_SPLICE_SITES(genome_gtf_ch)
        splice_sites_ch = EXTRACT_SPLICE_SITES.out.splice_sites
        exons_ch = EXTRACT_SPLICE_SITES.out.exons
    } else {
        // 创建空文件通道
        splice_sites_ch = Channel.fromPath("NO_FILE").map { file -> file("splice_sites.txt") }
        exons_ch = Channel.fromPath("NO_FILE").map { file -> file("exons.txt") }
    }
    
    // 构建索引
    BUILD_HISAT2_INDEX(genome_fasta_ch, splice_sites_ch, exons_ch)
}