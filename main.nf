#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 定义工作流参数
// --- 输入参数 ---
params.srr_ids = ""                          // 输入的 SRR ID 列表，用逗号分隔
params.local_genome_path = ""                // 本地基因组文件路径
params.local_gtf_path = ""                   // 本地GTF文件路径
params.download_genome_url = ""              // 参考基因组下载 URL
params.download_gtf_url = ""                 // GTF 注释文件下载 URL

// --- 本地 FASTQ 文件参数 ---
params.local_fastq_files = ""           // 本地FASTQ文件路径列表，支持通配符

// --- 输出目录参数 ---
params.data = "./data"                       // 唯一的基础路径参数

// --- 流程控制参数 (默认均为 false) ---
params.run_download_srr = false              // 是否执行 SRR 数据下载
params.run_download_genome = false           // 是否执行参考基因组下载
params.run_build_star_index = false          // 是否执行 STAR 索引构建
params.run_fastp = false                     // 是否执行 fastp 质控
params.run_star_align = false                // 是否执行 STAR 比对
params.run_featurecounts = false             // 是否执行 featureCounts 定量



// 定义日志处理器
log.info """
================================================
RNA-seq 分析工作流（参数化控制版）
================================================
SRR ID 列表: ${params.srr_ids ?: "未提供"}
本地基因组路径: ${params.local_genome_path ?: "未提供"}
本地GTF路径: ${params.local_gtf_path ?: "未提供"}
基因组下载URL: ${params.download_genome_url ?: "未提供"}
GTF下载URL: ${params.download_gtf_url ?: "未提供"}
本地FASTQ文件列表: ${params.local_fastq_files ?: "未提供"}
主输出目录: ${params.data}/results
------------------------------------------------
流程控制:
  SRR 数据下载: ${params.run_download_srr}
  参考基因组下载: ${params.run_download_genome}
  STAR 索引构建: ${params.run_build_star_index}
  Fastp 质控: ${params.run_fastp}
  STAR 比对: ${params.run_star_align}
  FeatureCounts 定量: ${params.run_featurecounts}
------------------------------------------------
================================================
"""

// 辅助函数定义
def extractSampleId(file) {
    // 从文件路径提取样本标识符
    return file.baseName.replaceAll(/[._](R?[12]|1|2)$/, '')
}

def cleanSampleId(sample_id) {
    // 清理样本标识符，确保符合文件命名规范
    return sample_id
        .replaceAll(/[^a-zA-Z0-9_-]/, '_')  // 替换特殊字符
        .replaceAll(/_{2,}/, '_')           // 合并多个下划线
        .replaceAll(/^_+|_+$/, '')         // 移除首尾下划线
        .isEmpty() ? "sample_${System.currentTimeMillis()}" : sample_id
}

// 流程定义

process download_srr_data {
    tag "下载 SRR 数据: ${srr_id}"
    label 'process_low'
    
    input:
    val srr_id
    
    output:
    path "${srr_id}_1.fastq.gz", emit: read1
    path "${srr_id}_2.fastq.gz", emit: read2, optional: true
    path "${srr_id}.single.fastq.gz", emit: read_single, optional: true
    path "${srr_id}.srr.downloaded", emit: status_file
    
    script:
    """
    # 创建临时目录
    mkdir -p ./sra_temp
    mkdir -p ./fastq_temp
    
    # 下载 SRA 数据
    conda run -n sra_env prefetch ${srr_id} -O ./sra_temp
    
    # 转换为 FASTQ 格式
    conda run -n sra_env fasterq-dump ./sra_temp/${srr_id}.sra -O ./fastq_temp --split-files --threads ${task.cpus}
    
    # 检查是否为双端测序并压缩
    if [ -f "./fastq_temp/${srr_id}_2.fastq" ]; then
        gzip ./fastq_temp/${srr_id}_1.fastq
        gzip ./fastq_temp/${srr_id}_2.fastq
        mv ./fastq_temp/${srr_id}_1.fastq.gz ${srr_id}_1.fastq.gz
        mv ./fastq_temp/${srr_id}_2.fastq.gz ${srr_id}_2.fastq.gz
    else
        gzip ./fastq_temp/${srr_id}.fastq
        mv ./fastq_temp/${srr_id}.fastq.gz ${srr_id}.single.fastq.gz
    fi
    
    # 清理临时文件
    rm -rf ./sra_temp ./fastq_temp
    
    # 创建状态文件
    echo "SRR数据下载完成: ${srr_id}" > ${srr_id}.srr.downloaded
    """
    
    publishDir "${params.data}/fastq", mode: 'copy', pattern: "*.fastq.gz"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.downloaded"
}

process download_reference_genome {
    tag "下载参考基因组"
    label 'process_low'
    
    output:
    path "genome.fa", emit: genome_fasta
    path "annotation.gtf", emit: genome_gtf
    path "genome.downloaded", emit: status_file
    
    script:
    """
    # 下载基因组
    wget ${params.download_genome_url} -O genome.fa.gz
    gunzip genome.fa.gz
    
    # 下载 GTF
    wget ${params.download_gtf_url} -O annotation.gtf.gz
    gunzip annotation.gtf.gz
    
    # 创建状态文件
    touch "genome.downloaded"
    """
    
    publishDir "${file(params.local_genome_path)}", mode: 'copy'
    publishDir "${file(params.local_genome_path)}/logs", mode: 'copy', pattern: "*.downloaded"
    
}

process build_star_index {
    tag "构建 STAR 索引"
    label 'process_medium'
    
    input:
    path genome_fasta
    path genome_gtf
    
    output:
    path "star_index", emit: index_dir
    path "star_index.built", emit: status_file
    
    script:
    index_output_dir = "star_index"
    """
    mkdir -p ${index_output_dir}
    
    # 构建STAR索引
    conda run -n align_env STAR \\
        --runMode genomeGenerate \\
        --genomeDir ${index_output_dir} \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${genome_gtf} \\
        --sjdbOverhang 100 \\
        --runThreadN ${task.cpus}
        
    # 创建状态文件
    touch "star_index.built"
    """
    
    publishDir "${file(params.local_genome_path).getParent()}", mode: 'copy', overwrite: false, pattern: "star_index"
    publishDir "${file(params.local_genome_path).getParent()}/logs", mode: 'copy', pattern: "*.built"
}

process run_fastp {
    tag "Fastp 质控: ${srr_id}"
    label 'process_medium'
    
    input:
    tuple val(srr_id), path(read1), path(read2), path(read_single)
    
    output:
    tuple val(srr_id), path("${srr_id}_1.trimmed.fastq.gz"), path("${srr_id}_2.trimmed.fastq.gz"), emit: trimmed_paired_reads, optional: true
    tuple val(srr_id), path("${srr_id}.single.trimmed.fastq.gz"), emit: trimmed_single_read, optional: true
    path "${srr_id}.fastp.html", emit: html_report
    path "${srr_id}.fastp.json", emit: json_report
    path "${srr_id}.fastp.done", emit: status_file
    
    script:
    """
    # 确定输入文件
    if [ -f "${read_single}" ]; then
        # 单端测序
        conda run -n qc_env fastp \\
            -i ${read_single} \\
            -o ${srr_id}.single.trimmed.fastq.gz \\
            --html ${srr_id}.fastp.html \\
            --json ${srr_id}.fastp.json \\
            --thread ${task.cpus}
    else
        # 双端测序
        conda run -n qc_env fastp \\
            -i ${read1} \\
            -I ${read2} \\
            -o ${srr_id}_1.trimmed.fastq.gz \\
            -O ${srr_id}_2.trimmed.fastq.gz \\
            --html ${srr_id}.fastp.html \\
            --json ${srr_id}.fastp.json \\
            --thread ${task.cpus}
    fi
    
    # 创建状态文件
    touch "${srr_id}.fastp.done"
    """
    
    publishDir "${params.data}/results/fastp/${srr_id}", mode: 'copy', pattern: "*.{html,json,fastq.gz}"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.fastp.done"
}

process run_star_align {
    tag "STAR 比对: ${srr_id}"
    label 'process_high'
    
    input:
    tuple val(srr_id), path(read1), path(read2), path(read_single)
    path index_dir
    
    output:
    tuple val(srr_id), path("${srr_id}.bam"), emit: bam
    path "${srr_id}.bam.done", emit: status_file
    
    script:
    """
    # 确定输入文件
    if [ -f "${read_single}" ]; then
        # 单端测序
        conda run -n align_env STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${index_dir} \\
            --readFilesIn ${read_single} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix . \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard
    else
        # 双端测序
        conda run -n align_env STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${index_dir} \\
            --readFilesIn ${read1} ${read2} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix . \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard
    fi
    
    # 重命名输出文件
    mv Aligned.sortedByCoord.out.bam ${srr_id}.bam
    
    # 创建BAM索引
    conda run -n align_env samtools index ${srr_id}.bam
    
    # 清理临时文件 (可选)
    rm -f Log.out Log.progress.out SJ.out.tab
    
    # 创建状态文件
    touch "${srr_id}.bam.done"
    """
    
    publishDir "${params.data}/results/bam/${srr_id}", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.bam.done"
}

// 新增：批量运行featureCounts的流程
process run_featurecounts_batch {
    tag "FeatureCounts 批量定量"
    label 'process_medium'
    
    input:
    path bam_files // 这是一个包含所有BAM文件路径的列表
    path gtf_file
    
    output:
    path "all_samples.counts.txt", emit: counts_matrix
    path "all_samples.counts.txt.summary", emit: counts_summary
    path "featurecounts.done", emit: status_file
    
    script:
    // 将bam_files列表转换为空格分隔的字符串
    bam_files_str = bam_files.join(' ')
    
    """
    # 运行featureCounts
    conda run -n quant_env featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf_file} \\
        -o all_samples.counts.txt \\
        ${bam_files_str}
        
    # 创建状态文件
    touch "featurecounts.done"
    """
    
    publishDir "${params.data}/results/featurecounts", mode: 'copy'
    publishDir "${params.data}/logs", mode: 'copy', pattern: "*.done"
}


// 主工作流 - 根据参数选择性执行步骤
workflow {
    // 创建SRR ID输入通道
    srr_ids_ch = params.run_download_srr ?
        Channel.fromList(params.srr_ids.split(',').collect { it.trim() }.unique()) :
        Channel.empty()
    // 1. 处理参考基因组 - 使用Channel.empty() + mix模式
    downloaded_genome_ch = params.run_download_genome ?
        download_reference_genome() :
        Channel.empty()
    
    local_genome_ch = params.run_download_genome ?
        Channel.empty() :
        Channel.of([
            genome_fasta: file(params.local_genome_path),
            genome_gtf: file(params.local_gtf_path)
        ])
    
    // 合并下载和本地基因组Channel
    genome_files_ch = downloaded_genome_ch.mix(local_genome_ch)

    // 2. 构建 STAR 索引 - 使用Channel.empty() + mix模式
    built_index_ch = params.run_build_star_index ?
        genome_files_ch.map { genome -> build_star_index(genome.genome_fasta, genome.genome_gtf) }.flatten() :
        Channel.empty()
    
    local_index_ch = params.run_build_star_index ?
        Channel.empty() :
        Channel.of([
            index_dir: file(params.local_genome_path ?
                file(params.local_genome_path).getParent() + "/star_index" :
                "./data/genomes/star_index")
        ])
    
    // 合并构建和本地索引Channel
    star_index_ch = built_index_ch.mix(local_index_ch)

    // 3. 准备输入数据 - 使用Channel.empty() + mix模式处理SRR下载
    // SRR下载数据Channel
    srr_fastq_ch = params.run_download_srr ?
        download_srr_data(srr_ids_ch) :
        Channel.empty()
    
    // 组织SRR下载数据为统一格式
    srr_organized_ch = params.run_download_srr ?
        srr_fastq_ch.read1.map { read1 ->
            [ read1.baseName.replace('_1.fastq', ''), read1, null, null ]
        }.mix(
            srr_fastq_ch.read2.map { read2 ->
                [ read2.baseName.replace('_2.fastq', ''), null, read2, null ]
            }
        ).mix(
            srr_fastq_ch.read_single.map { read_single ->
                [ read_single.baseName.replace('.single.fastq', ''), null, null, read_single ]
            }
        ).groupTuple().map { srr_id, r1, r2, rs ->
            [srr_id, r1[0], r2[0], rs[0]]
        } :
        Channel.empty()
    
    // 本地FASTQ文件Channel（统一使用 local_fastq_files 参数）
    local_fastq_ch = (params.local_fastq_files && !params.run_download_srr) ?
        Channel.fromPath(params.local_fastq_files.split(','))
            .filter { file -> file.exists() }
            .map { file ->
                def sample_id = extractSampleId(file)
                def cleaned_id = cleanSampleId(sample_id)
                // 判断是否为单端或双端测序文件
                if (file.name.contains('_1.') || file.name.contains('_R1')) {
                    return [cleaned_id, file, null, null] // 双端测序的 read1
                } else if (file.name.contains('_2.') || file.name.contains('_R2')) {
                    return [cleaned_id, null, file, null] // 双端测序的 read2
                } else {
                    return [cleaned_id, null, null, file] // 单端测序
                }
            }
            .groupTuple()
            .map { sample_id, r1_list, r2_list, single_list ->
                def r1 = r1_list.find { it != null }
                def r2 = r2_list.find { it != null }
                def single = single_list.find { it != null }
                return [sample_id, r1, r2, single]
            } :
        Channel.empty()
    
    // 合并所有输入数据源
    raw_reads_ch = srr_organized_ch.mix(local_fastq_ch)
    
    // 4. 运行 Fastp（如果启用）
    if (params.run_fastp) {
        trimmed_reads_ch = run_fastp(raw_reads_ch)
    } else {
        // 如果不运行fastp，直接使用原始数据
        trimmed_reads_ch = raw_reads_ch
    }

    // 5. 运行 STAR 比对 - 使用Channel.empty() + mix模式
    star_align_input_ch = params.run_star_align ?
        trimmed_reads_ch.combine(star_index_ch.map { it.index_dir }) :
        Channel.empty()
    
    bam_files_ch = params.run_star_align ?
        run_star_align(star_align_input_ch) :
        Channel.empty()

    // 6. 运行 FeatureCounts - 使用Channel.empty() + mix模式
    featurecounts_input_ch = params.run_featurecounts ?
        bam_files_ch.bam.collect()
            .combine(genome_files_ch.map { it.genome_gtf }) :
        Channel.empty()
    
    featurecounts_results_ch = params.run_featurecounts ?
        run_featurecounts_batch(featurecounts_input_ch) :
        Channel.empty()
}


// 工作流完成时的处理
workflow.onComplete {
    log.info """
    工作流完成！
    
    执行时间: ${workflow.duration}
    成功: ${workflow.success}
    主输出目录: ${params.data}/results
    日志目录: ${params.data}/logs
    """
}

// 错误处理
workflow.onError {
    log.error """
    工作流执行失败！
    
    错误信息: ${workflow.errorMessage}
    """
    
    // 将错误信息写入文件
    def errorFile = file("${params.data}/logs/error.txt")
    errorFile.getParent().mkdirs()
    errorFile.text = """
    工作流执行失败！
    
    错误信息: ${workflow.errorMessage}
    """
}