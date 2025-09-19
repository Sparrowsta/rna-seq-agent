#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// FastP 质控工具的独立流水线
// 支持动态参数配置和智能优化

// === FastP 参数配置 ===
params.sample_groups = []                    // 样本配对信息
params.data = "."                           // 数据目录
params.results_dir = "${params.data}/results" // 时间戳结果目录 - 由prepare_node动态设置

// FastP 质控参数 - 默认模板（可被Agent智能优化）
params.qualified_quality_phred = 20         // 质量阈值 (Q20) - 默认适中值
params.unqualified_percent_limit = 40       // 低质量碱基比例限制 (40%) - 标准值
params.n_base_limit = 5                     // N碱基数量限制 - 严格控制
params.length_required = 15                 // 最小read长度 - RNA-seq推荐值
params.adapter_trimming = true              // 默认启用adapter修剪
params.quality_filtering = true             // 默认启用质量过滤
params.length_filtering = true              // 默认启用长度过滤

// 高级参数 - 硬编码最佳实践
// 若未由Agent显式指定，以下多数参数维持空值或合理默认，脚本中按需拼接
params.phred64 = null                       // -6/--phred64 输入质量编码
params.reads_to_process = null              // --reads_to_process 限制处理reads数
params.fix_mgi_id = null                    // --fix_mgi_id 修复MGI测序ID

params.detect_adapter_for_pe = (params.detect_adapter_for_pe != null ? params.detect_adapter_for_pe : true)

// 前后定长修剪与最大长度
params.trim_front1 = null
params.trim_tail1 = null
params.max_len1   = null
params.trim_front2 = null
params.trim_tail2 = null
params.max_len2   = null

// polyG / polyX 修剪
params.trim_poly_g = (params.trim_poly_g != null ? params.trim_poly_g : true)
params.poly_g_min_len = null
params.disable_trim_poly_g = (params.disable_trim_poly_g != null ? params.disable_trim_poly_g : false)
params.trim_poly_x = (params.trim_poly_x != null ? params.trim_poly_x : false)
params.poly_x_min_len = null

// 滑窗切除开关（布尔值）
params.cut_front = (params.cut_front != null ? params.cut_front : false)  // -5 启用前端切除
params.cut_tail  = (params.cut_tail  != null ? params.cut_tail  : false)  // -3 启用后端切除  
params.cut_right = (params.cut_right != null ? params.cut_right : false)  // -r 启用右端切除
params.cut_window_size = (params.cut_window_size != null ? params.cut_window_size : 4)
params.cut_mean_quality = (params.cut_mean_quality != null ? params.cut_mean_quality : 20)
params.cut_front_window_size = null
params.cut_front_mean_quality = null
params.cut_tail_window_size = null
params.cut_tail_mean_quality = null
params.cut_right_window_size = null
params.cut_right_mean_quality = null

// 质量/长度过滤细化
params.average_qual = null
params.disable_length_filtering = (params.disable_length_filtering != null ? params.disable_length_filtering : false)
params.length_limit = null
params.low_complexity_filter = (params.low_complexity_filter != null ? params.low_complexity_filter : false)
params.complexity_threshold = null

// PE 重叠校正与检测
params.correction = (params.correction != null ? params.correction : false)
params.overlap_len_require = null
params.overlap_diff_limit = null
params.overlap_diff_percent_limit = null

// 过表达序列分析
params.overrepresentation_analysis = (params.overrepresentation_analysis != null ? params.overrepresentation_analysis : true)
params.overrepresentation_sampling = null

// FastP 质控process - 使用硬编码最佳实践模板
process fastp_quality_control {
    cpus (params.resources?.fastp?.cpus ?: 4)
    memory (params.resources?.fastp?.memory ?: '4 GB')
    tag "FastP质控: ${sample_id}"
    
    publishDir "${params.results_dir}/fastp/${sample_id}", mode: 'copy', pattern: "*.{html,json,fastq.gz}"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}*trimmed.fastq.gz"), emit: qc_reads, optional: true
    path "${sample_id}.fastp.*", emit: qc_reports, optional: true
    
    script:
    // 构建标准化fastp命令模板 - 硬编码最佳实践
    if (read1 && read1.name != "NO_FILE" && read2 && read2.name != "NO_FILE") {
        // 双端测序 - 完整参数模板
        """
        echo "开始FastP质控: ${sample_id} (双端测序)"
        echo "使用参数: Q${params.qualified_quality_phred}, 最小长度${params.length_required}bp"
        
        micromamba run -n qc_env fastp \\
            --in1 ${read1} \\
            --in2 ${read2} \\
            --out1 ${sample_id}_1.trimmed.fastq.gz \\
            --out2 ${sample_id}_2.trimmed.fastq.gz \\
            --html ${sample_id}.fastp.html \\
            --json ${sample_id}.fastp.json \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.qualified_quality_phred} \\
            --unqualified_percent_limit ${params.unqualified_percent_limit} \\
            --n_base_limit ${params.n_base_limit} \\
            --length_required ${params.length_required} \\
            ${params.phred64 ? '--phred64' : ''} \\
            ${params.reads_to_process != null ? "--reads_to_process ${params.reads_to_process}" : ''} \\
            ${params.fix_mgi_id ? '--fix_mgi_id' : ''} \\
            ${params.trim_front1 != null ? "--trim_front1 ${params.trim_front1}" : ''} \\
            ${params.trim_tail1  != null ? "--trim_tail1 ${params.trim_tail1}"   : ''} \\
            ${params.max_len1    != null ? "--max_len1 ${params.max_len1}"       : ''} \\
            ${params.trim_front2 != null ? "--trim_front2 ${params.trim_front2}" : ''} \\
            ${params.trim_tail2  != null ? "--trim_tail2 ${params.trim_tail2}"   : ''} \\
            ${params.max_len2    != null ? "--max_len2 ${params.max_len2}"       : ''} \\
            ${params.cut_front ? '--cut_front' : ''} \\
            ${params.cut_tail ? '--cut_tail' : ''} \\
            ${params.cut_right ? '--cut_right' : ''} \\
            ${params.cut_front_window_size != null ? "--cut_front_window_size ${params.cut_front_window_size}" : ''} \\
            ${params.cut_front_mean_quality != null ? "--cut_front_mean_quality ${params.cut_front_mean_quality}" : ''} \\
            ${params.cut_tail_window_size != null ? "--cut_tail_window_size ${params.cut_tail_window_size}" : ''} \\
            ${params.cut_tail_mean_quality != null ? "--cut_tail_mean_quality ${params.cut_tail_mean_quality}" : ''} \\
            ${params.cut_right_window_size != null ? "--cut_right_window_size ${params.cut_right_window_size}" : ''} \\
            ${params.cut_right_mean_quality != null ? "--cut_right_mean_quality ${params.cut_right_mean_quality}" : ''} \\
            --cut_window_size ${params.cut_window_size} \\
            --cut_mean_quality ${params.cut_mean_quality} \\
            ${params.adapter_trimming ? '' : '--disable_adapter_trimming'} \\
            ${params.quality_filtering ? '' : '--disable_quality_filtering'} \\
            ${params.length_filtering && !params.disable_length_filtering ? '' : '--disable_length_filtering'} \\
            ${params.average_qual != null ? "--average_qual ${params.average_qual}" : ''} \\
            ${params.length_limit != null ? "--length_limit ${params.length_limit}" : ''} \\
            ${params.low_complexity_filter ? '--low_complexity_filter' : ''} \\
            ${params.complexity_threshold != null ? "--complexity_threshold ${params.complexity_threshold}" : ''} \\
            ${params.disable_trim_poly_g ? '--disable_trim_poly_g' : (params.trim_poly_g ? '--trim_poly_g' : '')} \\
            ${params.poly_g_min_len != null ? "--poly_g_min_len ${params.poly_g_min_len}" : ''} \\
            ${params.trim_poly_x ? '--trim_poly_x' : ''} \\
            ${params.poly_x_min_len != null ? "--poly_x_min_len ${params.poly_x_min_len}" : ''} \\
            ${params.overrepresentation_analysis ? '--overrepresentation_analysis' : ''} \\
            ${params.overrepresentation_sampling != null ? "--overrepresentation_sampling ${params.overrepresentation_sampling}" : ''} \\
            ${params.detect_adapter_for_pe ? '--detect_adapter_for_pe' : ''} \\
            ${params.correction ? '--correction' : ''} \\
            ${params.overlap_len_require != null ? "--overlap_len_require ${params.overlap_len_require}" : ''} \\
            ${params.overlap_diff_limit != null ? "--overlap_diff_limit ${params.overlap_diff_limit}" : ''} \\
            ${params.overlap_diff_percent_limit != null ? "--overlap_diff_percent_limit ${params.overlap_diff_percent_limit}" : ''} \\
            --verbose
        
        echo "FastP质控完成: ${sample_id}"
        """
    } else if (read1 && read1.name != "NO_FILE") {
        // 单端测序 - 完整参数模板
        """
        echo "开始FastP质控: ${sample_id} (单端测序)"
        echo "使用参数: Q${params.qualified_quality_phred}, 最小长度${params.length_required}bp"
        
        micromamba run -n qc_env fastp \\
            --in1 ${read1} \\
            --out1 ${sample_id}.single.trimmed.fastq.gz \\
            --html ${sample_id}.fastp.html \\
            --json ${sample_id}.fastp.json \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.qualified_quality_phred} \\
            --unqualified_percent_limit ${params.unqualified_percent_limit} \\
            --n_base_limit ${params.n_base_limit} \\
            --length_required ${params.length_required} \\
            ${params.phred64 ? '--phred64' : ''} \\
            ${params.reads_to_process != null ? "--reads_to_process ${params.reads_to_process}" : ''} \\
            ${params.fix_mgi_id ? '--fix_mgi_id' : ''} \\
            ${params.trim_front1 != null ? "--trim_front1 ${params.trim_front1}" : ''} \\
            ${params.trim_tail1  != null ? "--trim_tail1 ${params.trim_tail1}"   : ''} \\
            ${params.max_len1    != null ? "--max_len1 ${params.max_len1}"       : ''} \\
            ${params.cut_front ? '--cut_front' : ''} \\
            ${params.cut_tail ? '--cut_tail' : ''} \\
            ${params.cut_right ? '--cut_right' : ''} \\
            ${params.cut_front_window_size != null ? "--cut_front_window_size ${params.cut_front_window_size}" : ''} \\
            ${params.cut_front_mean_quality != null ? "--cut_front_mean_quality ${params.cut_front_mean_quality}" : ''} \\
            ${params.cut_tail_window_size != null ? "--cut_tail_window_size ${params.cut_tail_window_size}" : ''} \\
            ${params.cut_tail_mean_quality != null ? "--cut_tail_mean_quality ${params.cut_tail_mean_quality}" : ''} \\
            ${params.cut_right_window_size != null ? "--cut_right_window_size ${params.cut_right_window_size}" : ''} \\
            ${params.cut_right_mean_quality != null ? "--cut_right_mean_quality ${params.cut_right_mean_quality}" : ''} \\
            --cut_window_size ${params.cut_window_size} \\
            --cut_mean_quality ${params.cut_mean_quality} \\
            ${params.adapter_trimming ? '' : '--disable_adapter_trimming'} \\
            ${params.quality_filtering ? '' : '--disable_quality_filtering'} \\
            ${params.length_filtering && !params.disable_length_filtering ? '' : '--disable_length_filtering'} \\
            ${params.average_qual != null ? "--average_qual ${params.average_qual}" : ''} \\
            ${params.length_limit != null ? "--length_limit ${params.length_limit}" : ''} \\
            ${params.low_complexity_filter ? '--low_complexity_filter' : ''} \\
            ${params.complexity_threshold != null ? "--complexity_threshold ${params.complexity_threshold}" : ''} \\
            ${params.disable_trim_poly_g ? '--disable_trim_poly_g' : (params.trim_poly_g ? '--trim_poly_g' : '')} \\
            ${params.poly_g_min_len != null ? "--poly_g_min_len ${params.poly_g_min_len}" : ''} \\
            ${params.trim_poly_x ? '--trim_poly_x' : ''} \\
            ${params.poly_x_min_len != null ? "--poly_x_min_len ${params.poly_x_min_len}" : ''} \\
            ${params.overrepresentation_analysis ? '--overrepresentation_analysis' : ''} \\
            ${params.overrepresentation_sampling != null ? "--overrepresentation_sampling ${params.overrepresentation_sampling}" : ''} \\
            --verbose
        
        echo "FastP质控完成: ${sample_id}"
        """
    } else {
        error "无效的FASTQ文件配置: sample_id=${sample_id}, read1=${read1?.name}, read2=${read2?.name}"
    }
}

// 主工作流
workflow {
    // 创建样本Channel
    def sample_groups = params.sample_groups ?: []
    
    if (sample_groups.isEmpty()) {
        error "错误: 未提供sample_groups参数，请通过-params-file传入样本信息"
    }
    
    fastq_samples_ch = Channel.from(sample_groups)
        .map { group ->
            def sample_id = group.sample_id
            def read1 = group.read1 ? file(group.read1) : file("NO_FILE")
            def read2 = group.read2 ? file(group.read2) : file("NO_FILE")
            
            return [sample_id, read1, read2]
        }
    
    // 执行FastP质控
    fastp_quality_control(fastq_samples_ch)
}
