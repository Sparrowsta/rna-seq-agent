#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * FeatureCounts基因定量流水线（DSL2）
 * 参考 fastp.nf / star.nf 的参数与结构，适配工具层调用
 * 输入：STAR比对后的BAM文件（JSON: [{sample_id, bam_file}, ...]）
 * 输出：每样本计数文件 + 汇总矩阵
 */

// ---------------- 参数定义（与工具层对齐） ----------------
// JSON字符串：[{"sample_id": "S1", "bam_file": "/path/to/S1.bam"}, ...]
params.input_bam_list = ""
params.gtf_file = ""
params.results_dir = ""
params.work_dir = ""

// FeatureCounts 参数（默认值可被 -params-file 覆盖）
params.feature_type = params.feature_type ?: "exon"
params.attribute_type = params.attribute_type ?: "gene_id"
params.strand_specificity = (params.strand_specificity != null ? params.strand_specificity as int : 0)
params.min_mapping_quality = (params.min_mapping_quality != null ? params.min_mapping_quality as int : 10)
params.count_reads_pairs = (params.count_reads_pairs != null ? params.count_reads_pairs as boolean : false)
params.count_multi_mapping_reads = (params.count_multi_mapping_reads != null ? params.count_multi_mapping_reads as boolean : false)
params.ignore_duplicates = (params.ignore_duplicates != null ? params.ignore_duplicates as boolean : false)
params.require_both_ends_mapped = (params.require_both_ends_mapped != null ? params.require_both_ends_mapped as boolean : false)
params.exclude_chimeric = (params.exclude_chimeric != null ? params.exclude_chimeric as boolean : false)

// ---------------- 输入验证 ----------------
if (!params.input_bam_list) error "Missing required parameter: input_bam_list"
if (!params.gtf_file) error "Missing required parameter: gtf_file"
if (!params.results_dir) error "Missing required parameter: results_dir"


// ---------------- 进程：一次性处理全部样本 ----------------
process FEATURECOUNTS_ALL {
    cpus (params.resources?.featurecounts?.cpus ?: 4)
    memory (params.resources?.featurecounts?.memory ?: '16 GB')
    publishDir "${params.results_dir}/featurecounts", mode: 'copy'

    input:
    // 接收唯一命名的 BAM 列表（来自 STAR：<sample_id>.Aligned.sortedByCoord.out.bam）
    path bam_files
    path gtf_file
    val sample_ids

    output:
    path "all_samples.featureCounts", emit: matrix_all
    path "all_samples.featureCounts.summary", emit: summary_all
    path "*/*.featureCounts", emit: per_sample_counts, optional: true
    path "*/*.featureCounts.summary", emit: per_sample_summaries, optional: true
    path "merged_counts_matrix.txt", emit: matrix_compat

    script:
    // 使用原始 BAM 文件名（已包含 sample_id 前缀），避免任何内部软链接
    def sample_count = sample_ids.size()
    def bams_arg = bam_files.collect { it.name }.join(' ')
    """
    set -e
    echo "一次性执行 featureCounts，样本数: ${sample_count}"

    micromamba run -n quant_env featureCounts \\
        -a ${gtf_file} \\
        -o all_samples.featureCounts \\
        -t ${params.feature_type} \\
        -g ${params.attribute_type} \\
        -Q ${params.min_mapping_quality} \\
        -T ${task.cpus} \\
        -s ${params.strand_specificity} \\
        ${params.count_reads_pairs ? '-p' : ''} \\
        ${params.count_multi_mapping_reads ? '-M' : ''} \\
        ${params.ignore_duplicates ? '--ignoreDup' : ''} \\
        ${params.require_both_ends_mapped ? '-B' : ''} \\
        ${params.exclude_chimeric ? '-C' : ''} \\
        ${bams_arg}

    if [ -f "all_samples.featureCounts" ]; then
        cp all_samples.featureCounts merged_counts_matrix.txt
        echo "FeatureCounts批量计数完成"
    else
        echo "错误: FeatureCounts执行失败"
        touch merged_counts_matrix.txt
    fi
    """
}

// ---------------- 主工作流 ----------------
workflow {
    def sample_data = new groovy.json.JsonSlurper().parseText(params.input_bam_list)
    if (!(sample_data instanceof List) || sample_data.isEmpty()) {
        error "input_bam_list 解析为空或格式错误"
    }

    def bam_list = sample_data.collect { file(it.bam_file) }
    def sample_ids = sample_data.collect { it.sample_id as String }

    bam_files_ch = Channel.value(bam_list)
    gtf_file_ch = Channel.value(file(params.gtf_file))
    sample_ids_ch = Channel.value(sample_ids)

    FEATURECOUNTS_ALL(bam_files_ch, gtf_file_ch, sample_ids_ch)
}
