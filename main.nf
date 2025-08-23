#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// =====================================================
// 重构后的RNA-seq分析工作流
// 支持工具选择和本地文件检查的统一架构
// =====================================================

// === 核心配置参数 ===
params.genome_version = "hg38"              // 基因组版本，所有路径从genomes.json获取

// === 数据输入控制 ===
params.local_fastq_files = ""               // 本地FASTQ文件路径列表
params.sample_groups = ""                   // Agent分析的样本配对信息(JSON字符串)

// === 基因组和索引控制 ===  
params.run_download_genome = false          // 是否下载基因组
params.run_build_star_index = false         // 是否构建STAR索引

// === 工具选择参数（可扩展设计）===
params.qc_tool = "fastp"                    // 质控工具: fastp, trimmomatic, none
params.align_tool = "star"                  // 比对工具: star, hisat2, none  
params.quant_tool = "featurecounts"         // 定量工具: featurecounts, htseq, none

// === 输出目录 ===
params.data = "."

// =====================================================
// 辅助函数定义
// =====================================================

// 获取完整的基因组路径信息
def getGenomePaths(genome_version) {
    def config_file = file("${projectDir}/config/genomes.json")
    def json_text = config_file.text
    def config = new groovy.json.JsonSlurper().parseText(json_text)
    def genome_config = config[genome_version]
    
    return [
        fasta: file(genome_config.fasta_path),
        gtf: file(genome_config.gtf_path),
        fasta_url: genome_config.fasta_url,
        gtf_url: genome_config.gtf_url,
        star_index: file(genome_config.fasta_path).getParent().resolve("star_index"),
        species: genome_config.species,
        version: genome_config.version
    ]
}

// 获取基因组目录路径
def getGenomeDir(genome_version) {
    def genome_paths = getGenomePaths(genome_version)
    return "genomes/${genome_paths.species}/${genome_version}"
}


// =====================================================
// 日志信息
// =====================================================

// 获取基因组配置信息用于日志显示
def global_genome_paths = getGenomePaths(params.genome_version)

log.info """
================================================
RNA-seq 分析工作流（重构版 - 工具选择模式）
================================================
基因组版本: ${params.genome_version} (${global_genome_paths.species})
基因组文件: ${global_genome_paths.fasta}
GTF文件: ${global_genome_paths.gtf}
STAR索引: ${global_genome_paths.star_index}

本地FASTQ文件: ${params.local_fastq_files ?: "未提供"}
主输出目录: ${params.data}/results
------------------------------------------------
工具选择:
  质控工具: ${params.qc_tool}
  比对工具: ${params.align_tool}
  定量工具: ${params.quant_tool}
------------------------------------------------
执行控制:
  参考基因组下载: ${params.run_download_genome}
  STAR 索引构建: ${params.run_build_star_index}
================================================
"""

// 统一的基因组管理process
process prepare_genome_data {
    tag "准备基因组数据: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "*.{fa,gtf}", enabled: params.run_download_genome
    publishDir "${getGenomeDir(params.genome_version)}/logs", mode: 'copy', pattern: "*.prepared"
    
    output:
    path "${params.genome_version}.fa", emit: genome_fasta
    path "*.gtf", emit: genome_gtf
    path "genome.prepared", emit: status_file
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def fasta_path = local_genome_paths.fasta.toString()
    def gtf_path = local_genome_paths.gtf.toString()
    def gtf_filename = file(gtf_path).getName()
    
    if (params.run_download_genome) {
        // 下载基因组数据
        """
        # 下载基因组
        wget ${local_genome_paths.fasta_url} -O ${params.genome_version}.fa.gz
        gunzip ${params.genome_version}.fa.gz
        
        # 下载 GTF
        wget ${local_genome_paths.gtf_url} -O ${gtf_filename}.gz
        gunzip ${gtf_filename}.gz
        
        echo "基因组数据下载完成: ${params.genome_version}" > genome.prepared
        """
    } else {
        // 使用本地基因组文件
        """
        # 检查本地文件是否存在
        if [ ! -f "${fasta_path}" ]; then
            echo "错误: 本地基因组文件不存在: ${fasta_path}"
            exit 1
        fi
        
        if [ ! -f "${gtf_path}" ]; then
            echo "错误: 本地GTF文件不存在: ${gtf_path}"
            exit 1
        fi
        
        # 创建软链接
        ln -s \$(realpath "${fasta_path}") ${params.genome_version}.fa
        ln -s \$(realpath "${gtf_path}") ${gtf_filename}
        
        echo "本地基因组数据准备完成: ${params.genome_version}" > genome.prepared
        """
    }
}

// 统一的STAR索引构建process
process prepare_star_index {
    tag "准备STAR索引: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "star_index", enabled: params.run_build_star_index
    publishDir "${getGenomeDir(params.genome_version)}/logs", mode: 'copy', pattern: "*.prepared"
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "star_index", emit: index_dir
    path "star_index.prepared", emit: status_file
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def star_index_path = local_genome_paths.star_index.toString()
    
    if (params.run_build_star_index) {
        // 构建STAR索引
        """
        mkdir -p star_index
        
        micromamba run -n align_env STAR \\
            --runMode genomeGenerate \\
            --genomeDir star_index \\
            --genomeFastaFiles ${genome_fasta} \\
            --sjdbGTFfile ${genome_gtf} \\
            --sjdbOverhang 100 \\
            --runThreadN ${task.cpus}
        
        echo "STAR索引构建完成: ${params.genome_version}" > star_index.prepared
        """
    } else {
        // 使用本地STAR索引
        """
        # 检查本地索引是否存在
        if [ ! -d "${star_index_path}" ]; then
            echo "错误: 本地STAR索引不存在: ${star_index_path}"
            exit 1
        fi
        
        # 创建软链接
        ln -s \$(realpath "${star_index_path}") star_index
        
        echo "本地STAR索引准备完成: ${params.genome_version}" > star_index.prepared
        """
    }
}

// 统一的质控process - 支持工具选择
process run_quality_control {
    tag "质控: ${params.qc_tool} - ${sample_id}"
    
    publishDir "${params.data}/results/${params.qc_tool}/${sample_id}", mode: 'copy', pattern: "*.{html,json,fastq.gz}"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.done"
    
    input:
    tuple val(sample_id), path(read1), path(read2), path(read_single)
    
    output:
    tuple val(sample_id), path("${sample_id}_*trimmed.fastq.gz"), emit: qc_reads, optional: true
    path "${sample_id}.${params.qc_tool}.*", emit: qc_reports, optional: true
    path "${sample_id}.${params.qc_tool}.done", emit: status_file
    
    when:
    params.qc_tool != "none"
    
    script:
    // 根据工具选择直接执行相应逻辑
    if (params.qc_tool == "fastp") {
        if (read_single && read_single.name != "NO_FILE") {
            // 单端测序
            """
            micromamba run -n qc_env fastp \\
                -i ${read_single} \\
                -o ${sample_id}.single.trimmed.fastq.gz \\
                --html ${sample_id}.fastp.html \\
                --json ${sample_id}.fastp.json \\
                --thread ${task.cpus}
            
            touch ${sample_id}.fastp.done
            """
        } else if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
            // 双端测序
            """
            micromamba run -n qc_env fastp \\
                -i ${read1} \\
                -I ${read2} \\
                -o ${sample_id}_1.trimmed.fastq.gz \\
                -O ${sample_id}_2.trimmed.fastq.gz \\
                --html ${sample_id}.fastp.html \\
                --json ${sample_id}.fastp.json \\
                --thread ${task.cpus}
            
            touch ${sample_id}.fastp.done
            """
        } else {
            error "无效的FASTQ文件配置: sample_id=${sample_id}, read1=${read1?.name}, read2=${read2?.name}, single=${read_single?.name}"
        }
    } else {
        error "不支持的质控工具: ${params.qc_tool}"
    }
}

// 统一的比对process - 支持工具选择
process run_alignment {
    tag "比对: ${params.align_tool} - ${sample_id}"
    
    publishDir "${params.data}/results/bam/${sample_id}", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.done"
    
    input:
    tuple val(sample_id), path(reads)
    path index_dir
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam_files
    path "${sample_id}.${params.align_tool}.*", emit: align_logs, optional: true
    path "${sample_id}.${params.align_tool}.done", emit: status_file
    
    when:
    params.align_tool != "none"
    
    script:
    // 根据工具选择直接执行相应逻辑
    if (params.align_tool == "star") {
        // 重新组织reads参数
        def read1 = reads.find { it.name.contains('_1.') || it.name.contains('_R1.') || it.name.endsWith('_1.fastq.gz') || it.name.contains('1.trimmed') }
        def read2 = reads.find { it.name.contains('_2.') || it.name.contains('_R2.') || it.name.endsWith('_2.fastq.gz') || it.name.contains('2.trimmed') }
        def read_single = reads.find { it.name.contains('.single.') }
        
        if (read_single && read_single.name != "NO_FILE") {
            // 单端测序
            """
            micromamba run -n align_env STAR \\
                --runThreadN ${task.cpus} \\
                --genomeDir ${index_dir} \\
                --readFilesIn ${read_single} \\
                --readFilesCommand zcat \\
                --outFileNamePrefix ${sample_id}. \\
                --outSAMtype BAM SortedByCoordinate \\
                --outSAMunmapped Within \\
                --outSAMattributes Standard
            
            mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
            micromamba run -n align_env samtools index ${sample_id}.bam
            
            touch ${sample_id}.star.done
            """
        } else if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
            // 双端测序
            """
            micromamba run -n align_env STAR \\
                --runThreadN ${task.cpus} \\
                --genomeDir ${index_dir} \\
                --readFilesIn ${read1} ${read2} \\
                --readFilesCommand zcat \\
                --outFileNamePrefix ${sample_id}. \\
                --outSAMtype BAM SortedByCoordinate \\
                --outSAMunmapped Within \\
                --outSAMattributes Standard \\
                --outSAMstrandField intronMotif \\
                --outFilterIntronMotifs RemoveNoncanonical
            
            mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
            micromamba run -n align_env samtools index ${sample_id}.bam
            
            touch ${sample_id}.star.done
            """
        } else {
            error "无效的FASTQ文件配置: sample_id=${sample_id}, read1=${read1?.name}, read2=${read2?.name}, single=${read_single?.name}"
        }
    } else {
        error "不支持的比对工具: ${params.align_tool}"
    }
}

// 统一的定量process - 支持工具选择
process run_quantification {
    tag "定量: ${params.quant_tool}"
    
    publishDir "${params.data}/results/${params.quant_tool}", mode: 'copy'
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.done"
    
    input:
    path bam_files
    path gtf_file
    
    output:
    path "all_samples.counts.txt", emit: counts_matrix, optional: true
    path "all_samples.counts.txt.summary", emit: counts_summary, optional: true
    path "${params.quant_tool}.done", emit: status_file
    
    when:
    params.quant_tool != "none"
    
    script:
    // 根据工具选择直接执行相应逻辑
    if (params.quant_tool == "featurecounts") {
        def bam_files_str = bam_files.join(' ')
        """
        # 检测BAM文件是否为双端测序
        first_bam=\$(echo ${bam_files_str} | cut -d' ' -f1)
        is_paired=\$(micromamba run -n align_env samtools view -c -f 1 "\$first_bam")
        
        if [ "\$is_paired" -gt 0 ]; then
            echo "检测到双端测序数据，使用双端模式"
            micromamba run -n quant_env featureCounts \\
                -T ${task.cpus} \\
                -p \\
                -B \\
                -C \\
                -a ${gtf_file} \\
                -o all_samples.counts.txt \\
                ${bam_files_str}
        else
            echo "检测到单端测序数据，使用单端模式"
            micromamba run -n quant_env featureCounts \\
                -T ${task.cpus} \\
                -a ${gtf_file} \\
                -o all_samples.counts.txt \\
                ${bam_files_str}
        fi
        
        touch featurecounts.done
        """
    } else {
        error "不支持的定量工具: ${params.quant_tool}"
    }
}

// =====================================================
// 主工作流 - 本地数据处理流程
// =====================================================

workflow {
    // --- 1. 数据输入 ---
    // 使用Agent分析的样本配对信息，避免复杂的文件名解析
    def sample_groups_json = params.sample_groups
    if (!sample_groups_json) {
        error "缺少样本配对信息：params.sample_groups未提供。请通过Agent分析生成配对信息。"
    }
    
    // 解析Agent提供的样本分组信息
    def sample_groups = new groovy.json.JsonSlurper().parseText(sample_groups_json)
    
    // 创建样本Channel，使用Agent的配对分析结果
    fastq_samples_ch = Channel.from(sample_groups)
        .map { group ->
            def sample_id = group.sample_id
            def read1 = group.read1 ? file(group.read1) : file("NO_FILE")
            def read2 = group.read2 ? file(group.read2) : file("NO_FILE") 
            def single = (!group.read2 && group.read1) ? file(group.read1) : file("NO_FILE")
            
            return [sample_id, read1, read2, single]
        }

    // --- 2. 准备基因组数据 ---
    genome_data_ch = prepare_genome_data()
    
    // --- 3. 准备STAR索引 ---
    star_index_ch = prepare_star_index(
        genome_data_ch.genome_fasta,
        genome_data_ch.genome_gtf
    )
    
    // --- 4. 质控处理 ---
    if (params.qc_tool != "none") {
        qc_results_ch = run_quality_control(fastq_samples_ch)
        processed_reads_ch = qc_results_ch.qc_reads
    } else {
        processed_reads_ch = fastq_samples_ch
    }
    
    // --- 6. 序列比对 ---
    if (params.align_tool != "none") {
        align_results_ch = run_alignment(
            processed_reads_ch,
            star_index_ch.index_dir
        )
        bam_files_ch = align_results_ch.bam_files
    } else {
        bam_files_ch = Channel.empty()
    }
    
    // --- 7. 基因定量 ---
    if (params.quant_tool != "none" && params.align_tool != "none") {
        quant_results_ch = run_quantification(
            bam_files_ch.map { sample_id, bam -> bam }.collect(),
            genome_data_ch.genome_gtf
        )
    }
}

// =====================================================
// 工作流完成处理
// =====================================================

workflow.onComplete {
    log.info """
    ================================================
    工作流完成！
    ================================================
    执行时间: ${workflow.duration}
    成功: ${workflow.success}
    基因组版本: ${params.genome_version}
    使用的工具:
      - 质控: ${params.qc_tool}
      - 比对: ${params.align_tool}  
      - 定量: ${params.quant_tool}
    主输出目录: ${params.data}/results
    日志目录: ${params.data}/logs
    ================================================
    """
}

// 错误处理
workflow.onError {
    log.error """
    ================================================
    工作流执行失败！
    ================================================
    错误信息: ${workflow.errorMessage}
    基因组版本: ${params.genome_version}
    使用的工具:
      - 质控: ${params.qc_tool}
      - 比对: ${params.align_tool}
      - 定量: ${params.quant_tool}
    ================================================
    """
    
    // 将错误信息写入文件
    def errorFile = file("${params.data}/logs/error.txt")
    errorFile.getParent().mkdirs()
    errorFile.text = """
    工作流执行失败！
    
    错误信息: ${workflow.errorMessage}
    执行时间: ${workflow.duration}
    基因组版本: ${params.genome_version}
    工具配置: QC=${params.qc_tool}, Align=${params.align_tool}, Quant=${params.quant_tool}
    """
}