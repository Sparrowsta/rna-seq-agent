#!/usr/bin/env nextflow


// 定义工作流参数
// --- 输入参数 ---
params.srr_ids = ""                          // 输入的 SRR ID 列表，用逗号分隔
params.local_genome_path = ""                // 本地基因组文件路径
params.local_gtf_path = ""                   // 本地GTF文件路径
params.download_genome_url = ""              // 参考基因组下载 URL (例如, ftp://.../genome.fa.gz)
params.download_gtf_url = ""                 // GTF 注释文件下载 URL (例如, ftp://.../annotation.gtf.gz)

// --- 输出目录参数 ---
params.data = "./data"                       // 唯一的基础路径参数

// --- 流程控制参数 (默认均为 false) ---
params.run_download_srr = false              // 是否执行 SRR 数据下载
params.run_download_genome = false           // 是否执行参考基因组下载
params.run_build_star_index = false          // 是否执行 STAR 索引构建
params.run_fastp = false                     // 是否执行 fastp 质控
params.run_star_align = false                // 是否执行 STAR 比对
params.run_featurecounts = false             // 是否执行 featureCounts 定量

// --- 其他参数 ---
params.resume = false                        // 是否从上次中断的地方继续
params.star_overhang = 100                   // STAR sjdbOverhang 值
params.star_threads = 4                      // STAR 运行线程数
params.fastp_threads = 4                     // fastp 运行线程数
params.featurecounts_threads = 4             // featureCounts 运行线程数



// 将 SRR ID 字符串转换为列表，并创建输入通道
srr_ids_ch = Channel.empty()
if (params.run_download_srr) {
    srr_list = params.srr_ids.split(',').collect { it.trim() }.unique()
    srr_ids_ch = Channel.fromList(srr_list)
}

// 定义日志处理器
log.info """
================================================
RNA-seq 分析工作流（参数化控制版）
================================================
SRR ID 列表: ${params.srr_ids ?: "未提供"}
参考基因组 URL: ${params.genome_url ?: "未提供"}
GTF URL: ${params.gtf_url ?: "未提供"}
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
恢复模式: ${params.resume}
================================================
"""

// 流程定义

process download_srr_data {
    tag "下载 SRR 数据: ${srr_id}"
    label 'process_low'
    
    input:
    val srr_id from srr_ids_ch
    
    output:
    path "${srr_id}_1.fastq.gz", emit: read1
    path "${srr_id}_2.fastq.gz", emit: read2, optional: true
    path "${srr_id}.single.fastq.gz", emit: read_single, optional: true
    path "${params.status_dir}/${srr_id}.srr.downloaded", emit: status_file
    
    script:
    // 创建输出目录
    sra_output_dir = file("${params.sra_dir}/${srr_id}")
    fastq_output_dir = file("${params.fastq_dir}/${srr_id}")
    status_output_dir = file(params.status_dir)
    
    """
    mkdir -p ${sra_output_dir}
    mkdir -p ${fastq_output_dir}
    mkdir -p ${status_output_dir}
    
    # 下载 SRA 数据
    conda run -n sra_env prefetch ${srr_id} -O ${sra_output_dir}
    
    # 转换为 FASTQ 格式并压缩
    conda run -n sra_env fasterq-dump ${sra_output_dir}/${srr_id}.sra -O ${fastq_output_dir} --split-files --threads ${task.cpus}
    
    # 检查是否为双端测序
    if [ -f "${fastq_output_dir}/${srr_id}_2.fastq" ]; then
        mv ${fastq_output_dir}/${srr_id}_1.fastq ${fastq_output_dir}/${srr_id}_1.fastq.gz
        mv ${fastq_output_dir}/${srr_id}_2.fastq ${fastq_output_dir}/${srr_id}_2.fastq.gz
    else
        mv ${fastq_output_dir}/${srr_id}.fastq ${fastq_output_dir}/${srr_id}.single.fastq.gz
    fi
    
    # 清理 SRA 文件以节省空间
    rm -rf ${sra_output_dir}
    
    # 创建状态文件
    echo "SRR数据下载完成: ${srr_id}" > ${params.status_dir}/${srr_id}.srr.downloaded
    """
    
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
    
    publishDir "${params.data}/genome", mode: 'copy'
    publishDir "${params.data}/status", mode: 'copy', pattern: "*.downloaded"
    
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
        --sjdbOverhang ${params.star_overhang} \\
        --runThreadN ${params.star_threads}
        
    # 创建状态文件
    touch "star_index.built"
    """
    
    publishDir "${file(params.local_genome_path).getParent()}", mode: 'copy', overwrite: false, pattern: "star_index"
    publishDir "${params.data}/status", mode: 'copy', pattern: "*.built"
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
            --thread ${params.fastp_threads}
    else
        # 双端测序
        conda run -n qc_env fastp \\
            -i ${read1} \\
            -I ${read2} \\
            -o ${srr_id}_1.trimmed.fastq.gz \\
            -O ${srr_id}_2.trimmed.fastq.gz \\
            --html ${srr_id}.fastp.html \\
            --json ${srr_id}.fastp.json \\
            --thread ${params.fastp_threads}
    fi
    
    # 创建状态文件
    touch "${srr_id}.fastp.done"
    """
    
    publishDir "${params.data}/results/fastp", mode: 'copy'
    publishDir "${params.data}/status", mode: 'copy', pattern: "*.fastp.done"
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
            --runThreadN ${params.star_threads} \\
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
            --runThreadN ${params.star_threads} \\
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
    
    publishDir "${params.data}/results/bam", mode: 'copy'
    publishDir "${params.data}/status", mode: 'copy', pattern: "*.bam.done"
}

process run_featurecounts {
    tag "FeatureCounts 定量"
    label 'process_medium'
    
    input:
    path index_dir // 用于推断GTF文件位置
    tuple val(srr_id), path(bam_file)
    
    output:
    path "all_samples.counts.txt", emit: counts_matrix
    path "all_samples.counts.txt.summary", emit: counts_summary
    path "${params.status_dir}/featurecounts.done", emit: status_file
    
    script:
    featurecounts_output_dir = file(params.featurecounts_dir)
    status_output_dir = file(params.status_dir)
    // 假设GTF文件与索引在同一目录，或由用户通过params.gtf_path明确指定
    // 这里我们尝试从索引目录中找到GTF文件
    gtf_file_path = file("${index_dir}/../annotation.gtf") // 相对于索引目录的GTF路径
    if (!gtf_file_path.exists() && params.gtf_url) {
        // 如果索引目录下没有，尝试从主基因组目录找
        gtf_file_path = file("${params.genome_dir}/annotation.gtf")
    }
    
    """
    mkdir -p ${featurecounts_output_dir}
    mkdir -p ${status_output_dir}
    cd ${featurecounts_output_dir}
    
    # 创建一个包含所有BAM文件路径的临时文件
    # echo "${bam_file}" > bam_files_list.txt
    
    # 运行featureCounts
    conda run -n quant_env featureCounts \\
        -T ${params.featurecounts_threads} \\
        -a ${gtf_file_path} \\
        -o all_samples.counts.txt \\
        ${bam_file} # 如果bam_file是单个文件路径
        
    # 创建状态文件
    echo "featureCounts基因定量完成" > ${params.status_dir}/featurecounts.done
    """
    
    publishDir "${params.featurecounts_dir}", mode: 'copy', overwrite: false
}

// 主工作流 - 根据参数选择性执行步骤
workflow {
    // 1. 处理参考基因组
    if (params.run_download_genome) {
        // 下载参考基因组
        downloaded_genome = download_reference_genome()
        // 将下载的基因组路径设置为后续流程的参数
        params.genome_fasta = downloaded_genome.out.genome_fasta
        params.gtf_path = downloaded_genome.out.genome_gtf
    } else {
        // 直接使用提供的路径
        downloaded_genome = [ genome_fasta: file(params.local_genome_path), genome_gtf: file(params.local_gtf_path) ]
    }

    // 2. 构建 STAR 索引
    if (params.run_build_star_index) {
        star_index = build_star_index(downloaded_genome.out.genome_fasta, downloaded_genome.out.genome_gtf)
    } else {
        // 如果不构建索引，则假定索引存在于本地基因组文件旁边
        def star_index_path = file(params.local_genome_path).getParent() + "/star_index"
        if (!file(star_index_path).exists()) error "STAR索引目录未找到: ${star_index_path}。请运行索引构建或提供正确路径。"
        // 创建虚拟通道
        star_index = [ index_dir: file(star_index_path) ]
    }

    // 3. 下载 SRR 数据
    trimmed_reads_ch = Channel.empty()
    if (params.run_download_srr) {
        raw_fastq_ch = download_srr_data(srr_ids_ch)
        
        // 4. 运行 Fastp
        if (params.run_fastp) {
            // 重新组织 raw_fastq_ch 的输出以匹配 run_fastp 的输入
            // raw_fastq_ch 输出: read1, read2, read_single
            // run_fastp 输入: tuple val(srr_id), path(read1), path(read2), path(read_single)
            // 需要从文件名中提取 srr_id
            raw_fastq_ch.out.read1.map { read1 -> [ read1.baseName.replace('_1.fastq', ''), read1, null, null ] }.mix(
                raw_fastq_ch.out.read2.map { read2 -> [ read2.baseName.replace('_2.fastq', ''), null, read2, null ] }
            ).mix(
                raw_fastq_ch.out.read_single.map { read_single -> [ read_single.baseName.replace('.single.fastq', ''), null, null, read_single ] }
            ).groupTuple().map { srr_id, r1, r2, rs -> [srr_id, r1[0], r2[0], rs[0]] }.set { fastp_input_ch }
            
            trimmed_reads_ch = run_fastp(fastp_input_ch)
        } else {
            // 如果不运行fastp，则使用原始fastq作为比对输入
            // 需要将原始fastq的输出组织成与trimmed_reads_ch相同的格式
            raw_fastq_ch.out.read1.map { read1 -> [ read1.baseName.replace('_1.fastq', ''), read1, null, null ] }.mix(
                raw_fastq_ch.out.read2.map { read2 -> [ read2.baseName.replace('_2.fastq', ''), null, read2, null ] }
            ).mix(
                raw_fastq_ch.out.read_single.map { read_single -> [ read_single.baseName.replace('.single.fastq', ''), null, null, read_single ] }
            ).groupTuple().map { srr_id, r1, r2, rs -> [srr_id, r1[0], r2[0], rs[0]] }.set { trimmed_reads_ch }
        }
    }

    // 5. 运行 STAR 比对
    bam_files_ch = Channel.empty()
    if (params.run_star_align) {
        // 检查是否有reads数据（经过fastp或原始）
        if (trimmed_reads_ch) {
            // star_index.out.index_dir 是一个路径
            // trimmed_reads_ch 的输出格式需要匹配 run_star_align 的输入
            // run_star_align 输入: tuple val(srr_id), path(read1), path(read2), path(read_single), path(index_dir)
            // trimmed_reads_ch 输出: tuple val(srr_id), path(trimmed_read1), path(trimmed_read2), path(trimmed_single_read) (可选)
            
            // 重新组织 trimmed_reads_ch 的输出
            star_align_input_ch = trimmed_reads_ch.out.trimmed_paired_reads.map { srr_id, r1, r2 -> [srr_id, r1, r2, null] }.mix(
                trimmed_reads_ch.out.trimmed_single_read.map { srr_id, rs -> [srr_id, null, null, rs] }
            ).combine(star_index.out.index_dir) // 组合索引路径
            
            bam_files_ch = run_star_align(star_align_input_ch)
        } else {
            log.warn "跳过STAR比对，因为没有可用的reads数据。请确保SRR下载已启用并成功。"
        }
    }

    // 6. 运行 FeatureCounts
    if (params.run_featurecounts) {
        if (bam_files_ch) {
            // featureCounts 需要所有BAM文件和一个GTF文件
            // 当前 run_featurecounts 设计为处理单个BAM文件，这是不正确的。
            // 我们需要收集所有BAM文件，然后一次性运行featureCounts。
            
            // 收集所有BAM文件路径
            all_bams = bam_files_ch.out.bam.collect()
            
            // 将GTF文件路径也传递进去
            // 假设GTF路径来自 downloaded_genome 或 params.gtf_path
            gtf_for_featurecounts = params.run_download_genome ? downloaded_genome.out.genome_gtf : file(params.gtf_path)
            
            // 修改 run_featurecounts 以接收BAM文件列表和GTF文件
            // 由于原始的 run_featurecounts 输入定义不同，我们需要调整或创建一个新流程
            // 为了快速实现，我们在这里直接调用featureCounts，而不是通过单独的process
            // 这违背了Nextflow的原则，但为了快速修改...
            // 更好的做法是重构 run_featurecounts
            
            // 临时解决方案：在workflow块中直接执行shell命令（不推荐，但为了快速演示）
            // 或者，我们创建一个简化的 process 来处理这个情况。
            // 这里我们选择创建一个新的 process `run_featurecounts_batch`
            
            run_featurecounts_batch(all_bams, gtf_for_featurecounts)
        } else {
            log.warn "跳过FeatureCounts，因为没有可用的BAM文件。请确保STAR比对已启用并成功。"
        }
    }
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
        -T ${params.featurecounts_threads} \\
        -a ${gtf_file} \\
        -o all_samples.counts.txt \\
        ${bam_files_str}
        
    # 创建状态文件
    touch "featurecounts.done"
    """
    
    publishDir "${params.data}/results/featurecounts", mode: 'copy'
    publishDir "${params.data}/status", mode: 'copy', pattern: "*.done"
}


// 工作流完成时的处理
workflow.onComplete {
    log.info """
    工作流完成！
    
    执行时间: ${workflow.duration}
    成功: ${workflow.success}
    主输出目录: ${params.data}/results
    状态目录: ${params.data}/status
    """
}

// 错误处理
workflow.onError {
    log.error """
    工作流执行失败！
    
    错误信息: ${workflow.errorMessage}
    """
    
    // 将错误信息写入文件
    file("${params.data}/status/error.txt").text = """
    工作流执行失败！
    
    错误信息: ${workflow.errorMessage}
    """
}