#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// =====================================================
// 重构后的RNA-seq分析工作流
// 支持工具选择和本地文件检查的统一架构
// =====================================================

// === 核心配置参数 ===
params.genome_version = "hg38"              // 基因组版本，所有路径从genomes.json获取

// === 数据输入控制 ===
params.srr_ids = ""                         // SRR ID列表，用逗号分隔
params.local_fastq_files = ""               // 本地FASTQ文件路径，支持通配符
params.run_download_srr = true             // 是否下载SRR数据

// === 基因组和索引控制 ===  
params.run_download_genome = false          // 是否下载基因组
params.run_build_star_index = false         // 是否构建STAR索引

// === 工具选择参数（可扩展设计）===
params.qc_tool = "fastp"                    // 质控工具: fastp, trimmomatic, none
params.align_tool = "star"                  // 比对工具: star, hisat2, none  
params.quant_tool = "featurecounts"         // 定量工具: featurecounts, htseq, none

// === 输出目录 ===
params.data = "./data"

// =====================================================
// 辅助函数定义
// =====================================================

// 读取基因组配置的辅助函数
def loadGenomeConfig(genome_version) {
    def config_file = file("${projectDir}/config/genomes.json")
    if (!config_file.exists()) {
        error "基因组配置文件不存在: ${config_file}"
    }
    
    def json_text = config_file.text
    def config = new groovy.json.JsonSlurper().parseText(json_text)
    
    if (!config.containsKey(genome_version)) {
        error "不支持的基因组版本: ${genome_version}. 支持的版本: ${config.keySet().join(', ')}"
    }
    
    return config[genome_version]
}

// 获取完整的基因组路径信息
def getGenomePaths(genome_version) {
    def genome_config = loadGenomeConfig(genome_version)
    
    return [
        fasta: file(genome_config.fasta),
        gtf: file(genome_config.gtf),
        fasta_url: genome_config.fasta_url,
        gtf_url: genome_config.gtf_url,
        star_index: file(genome_config.fasta).getParent().resolve("star_index"),
        species: genome_config.species,
        version: genome_config.version
    ]
}

// 获取基因组目录路径
def getGenomeDir(genome_version) {
    def genome_config = loadGenomeConfig(genome_version)
    return "data/genomes/${genome_config.species}/${genome_version}"
}

// 检查本地质控结果
def checkLocalQCResults(sample_id, qc_tool) {
    def result_patterns = [
        "fastp": "${params.data}/results/fastp/${sample_id}/${sample_id}.fastp.json",
        "trimmomatic": "${params.data}/results/trimmomatic/${sample_id}/${sample_id}.trimmomatic.log",
        "none": ""
    ]
    
    def result_file = file(result_patterns[qc_tool] ?: "")
    return [
        exists: result_file.exists(),
        path: result_file,
        tool: qc_tool
    ]
}

// 检查本地比对结果  
def checkLocalAlignResults(sample_id, align_tool) {
    def result_patterns = [
        "star": "${params.data}/results/bam/${sample_id}/${sample_id}.bam",
        "hisat2": "${params.data}/results/bam/${sample_id}/${sample_id}.bam",
        "none": ""
    ]
    
    def result_file = file(result_patterns[align_tool] ?: "")
    def log_patterns = [
        "star": "${params.data}/results/bam/${sample_id}/Log.final.out",
        "hisat2": "${params.data}/results/bam/${sample_id}/${sample_id}.hisat2.log"
    ]
    def log_file = file(log_patterns[align_tool] ?: "")
    
    return [
        exists: result_file.exists() && log_file.exists(),
        bam_path: result_file,
        log_path: log_file,
        tool: align_tool
    ]
}

// 检查本地定量结果
def checkLocalQuantResults(quant_tool) {
    def result_patterns = [
        "featurecounts": "${params.data}/results/featurecounts/all_samples.counts.txt.summary",
        "htseq": "${params.data}/results/htseq/all_samples.counts.txt",
        "none": ""
    ]
    
    def result_file = file(result_patterns[quant_tool] ?: "")
    return [
        exists: result_file.exists(),
        path: result_file,
        tool: quant_tool
    ]
}

// 提取样本ID的辅助函数
def extractSampleId(file_path) {
    def file_name = file_path.getName()
    
    // 移除常见的测序文件后缀
    def clean_name = file_name
        .replaceAll(/\.fastq\.gz$/, '')
        .replaceAll(/\.fq\.gz$/, '')
        .replaceAll(/\.fastq$/, '')
        .replaceAll(/\.fq$/, '')
    
    // 移除双端测序标识
    clean_name = clean_name
        .replaceAll(/_[12]$/, '')
        .replaceAll(/_R[12]$/, '')
        .replaceAll(/\.single$/, '')
    
    return clean_name
}

// 清理样本ID
def cleanSampleId(sample_id) {
    return sample_id
        .replaceAll(/[^a-zA-Z0-9_-]/, '_')  // 替换特殊字符为下划线
        .replaceAll(/_+/, '_')              // 合并多个下划线
        .replaceAll(/^_|_$/, '')            // 移除首尾下划线
}

// Fastp脚本生成
def fastpScript(sample_id, read1, read2, read_single) {
    if (read_single && read_single.name != "NO_FILE") {
        // 单端测序
        return """
        conda run -n qc_env fastp \\
            -i ${read_single} \\
            -o ${sample_id}.single.trimmed.fastq.gz \\
            --html ${sample_id}.fastp.html \\
            --json ${sample_id}.fastp.json \\
            --thread ${task.cpus}
        
        touch ${sample_id}.fastp.done
        """
    } else if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
        // 双端测序
        return """
        conda run -n qc_env fastp \\
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
}

// STAR比对脚本生成
def starAlignScript(sample_id, read1, read2, read_single, index_dir) {
    if (read_single && read_single.name != "NO_FILE") {
        // 单端测序
        return """
        conda run -n align_env STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${index_dir} \\
            --readFilesIn ${read_single} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${sample_id}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard
        
        mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
        conda run -n align_env samtools index ${sample_id}.bam
        
        touch ${sample_id}.star.done
        """
    } else if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
        // 双端测序
        return """
        conda run -n align_env STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${index_dir} \\
            --readFilesIn ${read1} ${read2} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${sample_id}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard
        
        mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
        conda run -n align_env samtools index ${sample_id}.bam
        
        touch ${sample_id}.star.done
        """
    } else {
        error "无效的FASTQ文件配置: sample_id=${sample_id}, read1=${read1?.name}, read2=${read2?.name}, single=${read_single?.name}"
    }
}

// FeatureCounts脚本生成
def featureCountsScript(bam_files, gtf_file) {
    def bam_files_str = bam_files.join(' ')
    return """
    conda run -n quant_env featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf_file} \\
        -o all_samples.counts.txt \\
        ${bam_files_str}
    
    touch featurecounts.done
    """
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

SRR ID 列表: ${params.srr_ids ?: "未提供"}
本地FASTQ文件: ${params.local_fastq_files ?: "未提供"}
主输出目录: ${params.data}/results
------------------------------------------------
工具选择:
  质控工具: ${params.qc_tool}
  比对工具: ${params.align_tool}
  定量工具: ${params.quant_tool}
------------------------------------------------
执行控制:
  SRR 数据下载: ${params.run_download_srr}
  参考基因组下载: ${params.run_download_genome}
  STAR 索引构建: ${params.run_build_star_index}
================================================
"""

// =====================================================
// Process 定义
// =====================================================

// 统一的数据输入process - 整合SRR下载和本地FASTQ
process prepare_input_data {
    tag "准备输入数据"
    label 'process_low'
    
    publishDir "${params.data}/fastq", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.prepared"
    
    output:
    path "*.fastq.gz", emit: fastq_files, optional: true
    path "*.input.prepared", emit: status_files
    
    script:
    if (params.run_download_srr && params.srr_ids) {
        // 下载SRR数据
        def srr_list = params.srr_ids.split(',').collect { it.trim() }
        def srr_commands = srr_list.collect { srr_id ->
            """
            # 下载和处理 ${srr_id}
            mkdir -p ./sra_temp_${srr_id} ./fastq_temp_${srr_id}
            
            conda run -n sra_env prefetch ${srr_id} -O ./sra_temp_${srr_id}
            conda run -n sra_env fasterq-dump ./sra_temp_${srr_id}/${srr_id}/${srr_id}.sra -O ./fastq_temp_${srr_id} --split-files --threads ${task.cpus}
            
            if [ -f "./fastq_temp_${srr_id}/${srr_id}_2.fastq" ]; then
                gzip ./fastq_temp_${srr_id}/${srr_id}_1.fastq
                gzip ./fastq_temp_${srr_id}/${srr_id}_2.fastq
                mv ./fastq_temp_${srr_id}/${srr_id}_1.fastq.gz ${srr_id}_1.fastq.gz
                mv ./fastq_temp_${srr_id}/${srr_id}_2.fastq.gz ${srr_id}_2.fastq.gz
            else
                gzip ./fastq_temp_${srr_id}/${srr_id}.fastq
                mv ./fastq_temp_${srr_id}/${srr_id}.fastq.gz ${srr_id}.single.fastq.gz
            fi
            
            rm -rf ./sra_temp_${srr_id} ./fastq_temp_${srr_id}
            echo "SRR数据准备完成: ${srr_id}" > ${srr_id}.input.prepared
            """
        }.join('\n\n')
        
        return srr_commands
    } else if (params.local_fastq_files) {
        // 使用本地FASTQ文件
        """
        # 创建本地文件的软链接并重命名为标准格式
        for pattern in ${params.local_fastq_files.split(',').join(' ')}; do
            for file in ${projectDir}/\$pattern; do
                if [ -f "\$file" ]; then
                    # 获取文件名并创建软链接到当前工作目录
                    filename=\$(basename "\$file")
                    ln -s \$(realpath "\$file") "\$filename"
                    echo "链接文件: \$file -> \$filename"
                fi
            done
        done
        
        echo "本地FASTQ文件准备完成" > local_fastq.input.prepared
        """
    } else {
        error "必须提供SRR ID列表或本地FASTQ文件路径"
    }
}

// 统一的基因组管理process
process prepare_genome_data {
    tag "准备基因组数据: ${params.genome_version}"
    label 'process_low'
    
    publishDir "${getGenomeDir(params.genome_version)}", mode: 'copy', pattern: "*.{fa,gtf}", enabled: params.run_download_genome
    publishDir "${getGenomeDir(params.genome_version)}/logs", mode: 'copy', pattern: "*.prepared"
    
    output:
    path "genome.fa", emit: genome_fasta
    path "annotation.gtf", emit: genome_gtf
    path "genome.prepared", emit: status_file
    
    script:
    def local_genome_paths = getGenomePaths(params.genome_version)
    def fasta_path = local_genome_paths.fasta.toString()
    def gtf_path = local_genome_paths.gtf.toString()
    
    if (params.run_download_genome) {
        // 下载基因组数据
        """
        # 下载基因组
        wget ${local_genome_paths.fasta_url} -O genome.fa.gz
        gunzip genome.fa.gz
        
        # 下载 GTF
        wget ${local_genome_paths.gtf_url} -O annotation.gtf.gz
        gunzip annotation.gtf.gz
        
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
        ln -s \$(realpath "${fasta_path}") genome.fa
        ln -s \$(realpath "${gtf_path}") annotation.gtf
        
        echo "本地基因组数据准备完成: ${params.genome_version}" > genome.prepared
        """
    }
}

// 统一的STAR索引构建process
process prepare_star_index {
    tag "准备STAR索引: ${params.genome_version}"
    label 'process_medium'
    
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
        
        conda run -n align_env STAR \\
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
    label 'process_medium'
    
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
    // 检查本地结果文件
    def local_result = checkLocalQCResults(sample_id, params.qc_tool)
    
    if (local_result.exists) {
        // 使用本地文件
        """
        echo "使用本地质控结果: ${local_result.path}"
        ln -s ${local_result.path} .
        touch ${sample_id}.${params.qc_tool}.done
        """
    } else {
        // 根据工具选择执行相应逻辑
        if (params.qc_tool == "fastp") {
            if (read_single && read_single.name != "NO_FILE") {
                // 单端测序
                """
                conda run -n qc_env fastp \\
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
                conda run -n qc_env fastp \\
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
}

// 统一的比对process - 支持工具选择
process run_alignment {
    tag "比对: ${params.align_tool} - ${sample_id}"
    label 'process_high'
    
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
    // 检查本地结果文件
    def local_result = checkLocalAlignResults(sample_id, params.align_tool)
    
    if (local_result.exists) {
        // 使用本地文件
        """
        echo "使用本地比对结果: ${local_result.bam_path}"
        ln -s ${local_result.bam_path} ${sample_id}.bam
        ln -s ${local_result.log_path} .
        touch ${sample_id}.${params.align_tool}.done
        """
    } else {
        // 根据工具选择执行相应逻辑
        if (params.align_tool == "star") {
            // 需要重新组织reads参数
            def read1 = reads.find { it.name.contains('_1.') || it.name.contains('_R1') }
            def read2 = reads.find { it.name.contains('_2.') || it.name.contains('_R2') }
            def read_single = reads.find { it.name.contains('.single.') }
            
            if (read_single && read_single.name != "NO_FILE") {
                // 单端测序
                """
                conda run -n align_env STAR \\
                    --runThreadN ${task.cpus} \\
                    --genomeDir ${index_dir} \\
                    --readFilesIn ${read_single} \\
                    --readFilesCommand zcat \\
                    --outFileNamePrefix ${sample_id}. \\
                    --outSAMtype BAM SortedByCoordinate \\
                    --outSAMunmapped Within \\
                    --outSAMattributes Standard
                
                mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
                conda run -n align_env samtools index ${sample_id}.bam
                
                touch ${sample_id}.star.done
                """
            } else if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
                // 双端测序
                """
                conda run -n align_env STAR \\
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
                conda run -n align_env samtools index ${sample_id}.bam
                
                touch ${sample_id}.star.done
                """
            } else {
                error "无效的FASTQ文件配置: sample_id=${sample_id}, read1=${read1?.name}, read2=${read2?.name}, single=${read_single?.name}"
            }
        } else {
            error "不支持的比对工具: ${params.align_tool}"
        }
    }
}

// 统一的定量process - 支持工具选择
process run_quantification {
    tag "定量: ${params.quant_tool}"
    label 'process_medium'
    
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
    // 检查本地结果文件
    def local_result = checkLocalQuantResults(params.quant_tool)
    
    if (local_result.exists) {
        // 使用本地文件
        """
        echo "使用本地定量结果: ${local_result.path}"
        ln -s ${local_result.path} .
        touch ${params.quant_tool}.done
        """
    } else {
        // 根据工具选择执行相应逻辑
        if (params.quant_tool == "featurecounts") {
            def bam_files_str = bam_files.join(' ')
            """
            # 检测BAM文件是否为双端测序
            first_bam=\$(echo ${bam_files_str} | cut -d' ' -f1)
            is_paired=\$(conda run -n align_env samtools view -c -f 1 "\$first_bam")
            
            if [ "\$is_paired" -gt 0 ]; then
                echo "检测到双端测序数据，使用双端模式"
                conda run -n quant_env featureCounts \\
                    -T ${task.cpus} \\
                    -p \\
                    -B \\
                    -C \\
                    -a ${gtf_file} \\
                    -o all_samples.counts.txt \\
                    ${bam_files_str}
            else
                echo "检测到单端测序数据，使用单端模式"
                conda run -n quant_env featureCounts \\
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
}

// =====================================================
// 主工作流 - 简化的线性流程
// =====================================================

workflow {
    // 1. 准备输入数据（SRR下载或本地FASTQ）
    input_data_ch = prepare_input_data()
    
    // 2. 准备基因组数据（下载或本地）
    genome_data_ch = prepare_genome_data()
    
    // 3. 准备STAR索引（构建或本地）
    star_index_ch = prepare_star_index(
        genome_data_ch.genome_fasta,
        genome_data_ch.genome_gtf
    )
    
    // 4. 处理FASTQ文件，组织为样本格式
    fastq_samples_ch = input_data_ch.fastq_files
        .view { "Raw FASTQ files from input: $it" }
        .flatten()
        .view { "Flattened FASTQ file: $it" }
        .filter { it.name.endsWith('.fastq.gz') }
        .view { "Filtered FASTQ file: $it" }
        .map { file ->
            def sample_id = extractSampleId(file)
            def cleaned_id = cleanSampleId(sample_id)
            
            // 判断文件类型并返回样本信息
            if (file.name.contains('_1.') || file.name.contains('_R1')) {
                return [cleaned_id, 'read1', file]
            } else if (file.name.contains('_2.') || file.name.contains('_R2')) {
                return [cleaned_id, 'read2', file]
            } else {
                return [cleaned_id, 'single', file]
            }
        }
        .view { "Mapped FASTQ sample: $it" }
        .groupTuple()
        .map { sample_id, type_list, file_list ->
            def r1 = null
            def r2 = null
            def single = null
            
            // 根据类型分配文件
            for (int i = 0; i < type_list.size(); i++) {
                switch(type_list[i]) {
                    case 'read1':
                        r1 = file_list[i]
                        break
                    case 'read2':
                        r2 = file_list[i]
                        break
                    case 'single':
                        single = file_list[i]
                        break
                }
            }
            
            // 创建占位符文件对象以避免null值
            def null_file = file('NO_FILE')
            return [sample_id, r1 ?: null_file, r2 ?: null_file, single ?: null_file]
        }
        .view { "Final FASTQ samples for processing: $it" }
    
    // 5. 质控处理（如果启用）
    if (params.qc_tool != "none") {
        qc_results_ch = run_quality_control(fastq_samples_ch)
        processed_reads_ch = qc_results_ch.qc_reads
    } else {
        processed_reads_ch = fastq_samples_ch
    }
    
    // 6. 序列比对（如果启用）
    if (params.align_tool != "none") {
        align_results_ch = run_alignment(
            processed_reads_ch,
            star_index_ch.index_dir
        )
        bam_files_ch = align_results_ch.bam_files
    } else {
        bam_files_ch = Channel.empty()
    }
    
    // 7. 基因定量（如果启用）
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