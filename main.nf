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

// 进程级并行: 分离FASTA下载
process download_genome_fasta {
    tag "下载FASTA: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "*.fa", enabled: params.run_download_genome
    
    output:
    path "${params.genome_version}.fa", emit: genome_fasta
    
    when:
    params.run_download_genome
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    """
    echo "开始下载基因组FASTA: ${params.genome_version}" 
    wget ${local_genome_paths.fasta_url} -O ${params.genome_version}.fa.gz
    gunzip ${params.genome_version}.fa.gz
    echo "FASTA下载完成: ${params.genome_version}"
    """
}

// 进程级并行: 分离GTF下载
process download_genome_gtf {
    tag "下载GTF: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "*.gtf", enabled: params.run_download_genome
    
    output:
    path "${params.genome_version}.gtf", emit: genome_gtf
    
    when:
    params.run_download_genome
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def gtf_path = local_genome_paths.gtf.toString()
    def gtf_filename = file(gtf_path).getName()
    
    """
    echo "开始下载基因组GTF: ${params.genome_version}"
    wget ${local_genome_paths.gtf_url} -O ${params.genome_version}.gtf.gz
    gunzip ${params.genome_version}.gtf.gz
    echo "GTF下载完成: ${params.genome_version}"
    """
}

// 本地基因组文件准备
process prepare_local_genome {
    tag "准备本地基因组: ${params.genome_version}"
    
    output:
    path "${params.genome_version}.fa", emit: genome_fasta
    path "${params.genome_version}.gtf", emit: genome_gtf
    path "genome.prepared", emit: status_file
    
    when:
    !params.run_download_genome
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def fasta_path = local_genome_paths.fasta.toString()
    def gtf_path = local_genome_paths.gtf.toString()
    def gtf_filename = file(gtf_path).getName()
    
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
    ln -s \$(realpath "${gtf_path}") ${params.genome_version}.gtf
    
    echo "本地基因组数据准备完成: ${params.genome_version}" > genome.prepared
    """
}

// STAR索引构建process
process build_star_index {
    tag "构建STAR索引: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "star_index"
    publishDir "${getGenomeDir(params.genome_version)}/logs", mode: 'copy', pattern: "*.prepared"
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "star_index", emit: index_dir
    path "star_index.prepared", emit: status_file
    
    when:
    params.run_build_star_index
    
    script:
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
}

// STAR索引链接process
process link_star_index {
    tag "链接本地STAR索引: ${params.genome_version}"
    
    publishDir "${getGenomeDir(params.genome_version)}/logs", mode: 'copy', pattern: "*.prepared"
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "star_index", emit: index_dir
    path "star_index.prepared", emit: status_file
    
    when:
    !params.run_build_star_index
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def star_index_path = local_genome_paths.star_index.toString()
    
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
    // 解析Agent提供的样本分组信息
    def sample_groups = new groovy.json.JsonSlurper().parseText(params.sample_groups)
    
    // 创建样本Channel，使用Agent的配对分析结果
    fastq_samples_ch = Channel.from(sample_groups)
        .map { group ->
            def sample_id = group.sample_id
            def read1 = group.read1 ? file(group.read1) : file("NO_FILE")
            def read2 = group.read2 ? file(group.read2) : file("NO_FILE") 
            def single = (!group.read2 && group.read1) ? file(group.read1) : file("NO_FILE")
            
            return [sample_id, read1, read2, single]
        }

    // --- 2. 准备基因组数据 (when条件自动处理) ---
    download_genome_fasta()  // when: params.run_download_genome
    download_genome_gtf()    // when: params.run_download_genome  
    prepare_local_genome()   // when: !params.run_download_genome
    
    // 混合输出 - Nextflow自动处理空channel
    fasta_ch = download_genome_fasta.out.genome_fasta.mix(
        prepare_local_genome.out.genome_fasta
    )
    gtf_ch = download_genome_gtf.out.genome_gtf.mix(
        prepare_local_genome.out.genome_gtf
    )
    
    // --- 3. 准备STAR索引 (when条件自动处理) ---
    build_star_index(fasta_ch, gtf_ch)      // when: params.run_build_star_index
    link_star_index(fasta_ch, gtf_ch)       // when: !params.run_build_star_index
    
    // 混合两个索引process的输出
    star_index_ch = build_star_index.out.index_dir.mix(
        link_star_index.out.index_dir
    )
    
    // --- 4. 质控处理 (when条件自动处理) ---
    qc_results_ch = run_quality_control(fastq_samples_ch)  // when: params.qc_tool != "none"
    
    // 混合原始和处理后的reads
    processed_reads_ch = qc_results_ch.qc_reads.mix(
        fastq_samples_ch.filter { params.qc_tool == "none" }
    )
    
    // --- 5. 序列比对 (when条件自动处理) ---
    align_results_ch = run_alignment(
        processed_reads_ch,
        star_index_ch
    )  // when: params.align_tool != "none"
    
    // --- 6. 基因定量 (when条件自动处理) ---
    run_quantification(
        align_results_ch.bam_files.map { sample_id, bam -> bam }.collect(),
        gtf_ch
    )  // when: params.quant_tool != "none"
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