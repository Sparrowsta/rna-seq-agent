nextflow.enable.dsl=2

/*
 * STAR索引构建流水线
 * 输入：基因组FASTA和GTF文件
 * 输出：STAR索引目录
 * 支持所有STAR索引构建参数的动态传递
 */

// 基础必需参数
params.genome_fasta = ""
params.genome_gtf = ""
params.star_index_dir = ""

// STAR 索引构建参数 - 支持动态传递
params.sjdb_overhang = 100
params.sjdbOverhang = 100  // 兼容两种命名
params.runThreadN = 8
params.limitGenomeGenerateRAM = 32000000000

// 高级索引参数 - 可选
params.genomeSAindexNbases = null
params.genomeChrBinNbits = null
params.genomeSAsparseD = null
params.sjdbGTFfeatureExon = null
params.sjdbGTFtagExonParentTranscript = null
params.sjdbGTFtagExonParentGene = null
params.sjdbInsertSave = null

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

    # 构建STAR索引命令 - 支持动态参数
    STAR_CMD="micromamba run -n align_env STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles \$FASTA_FILE --sjdbGTFfile \$GTF_FILE"
    
    # 添加基础参数
    STAR_CMD="\$STAR_CMD --sjdbOverhang ${params.sjdbOverhang ?: params.sjdb_overhang}"
    STAR_CMD="\$STAR_CMD --runThreadN ${task.cpus}"
    STAR_CMD="\$STAR_CMD --limitGenomeGenerateRAM ${params.limitGenomeGenerateRAM}"
    
    # 添加可选高级参数（仅在非null时添加）
    if [[ "${params.genomeSAindexNbases}" != "null" && "${params.genomeSAindexNbases}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --genomeSAindexNbases ${params.genomeSAindexNbases}"
    fi
    
    if [[ "${params.genomeChrBinNbits}" != "null" && "${params.genomeChrBinNbits}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --genomeChrBinNbits ${params.genomeChrBinNbits}"
    fi
    
    if [[ "${params.genomeSAsparseD}" != "null" && "${params.genomeSAsparseD}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --genomeSAsparseD ${params.genomeSAsparseD}"
    fi
    
    if [[ "${params.sjdbGTFfeatureExon}" != "null" && "${params.sjdbGTFfeatureExon}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --sjdbGTFfeatureExon ${params.sjdbGTFfeatureExon}"
    fi
    
    if [[ "${params.sjdbGTFtagExonParentTranscript}" != "null" && "${params.sjdbGTFtagExonParentTranscript}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --sjdbGTFtagExonParentTranscript ${params.sjdbGTFtagExonParentTranscript}"
    fi
    
    if [[ "${params.sjdbGTFtagExonParentGene}" != "null" && "${params.sjdbGTFtagExonParentGene}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --sjdbGTFtagExonParentGene ${params.sjdbGTFtagExonParentGene}"
    fi
    
    if [[ "${params.sjdbInsertSave}" != "null" && "${params.sjdbInsertSave}" != "" ]]; then
        STAR_CMD="\$STAR_CMD --sjdbInsertSave ${params.sjdbInsertSave}"
    fi
    
    # 执行STAR索引构建
    echo "执行STAR索引构建命令: \$STAR_CMD"
    eval \$STAR_CMD
    
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