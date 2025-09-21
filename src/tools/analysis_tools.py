"""
RNA-seq智能分析助手 - 分析和报告工具

包含：
- parse_pipeline_results: 解析流水线结果数据并进行样本对齐整合
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.analysis")


@tool
def parse_pipeline_results(results_directory: str) -> Dict[str, Any]:
    """解析RNA-seq流水线的所有步骤结果并进行样本对齐整合

    Args:
        results_directory: 结果目录路径

    Returns:
        Dict包含所有步骤的解析结果、样本对齐数据和质量评估
    """
    try:
        from . import parse_fastp_results, parse_star_metrics, parse_hisat2_metrics, parse_featurecounts_metrics

        results_path = Path(results_directory)
        if not results_path.exists():
            return {
                "success": False,
                "error": f"结果目录不存在: {results_directory}"
            }

        logger.info(f"开始解析流水线结果: {results_directory}")

        # 解析各步骤结果
        fastp_result = parse_fastp_results.invoke({"results_directory": results_directory})
        star_result = parse_star_metrics.invoke({"results_directory": results_directory})
        hisat2_result = parse_hisat2_metrics.invoke({"results_directory": results_directory})
        fc_result = parse_featurecounts_metrics.invoke({"results_directory": results_directory})

        # 统计解析状态
        parsing_status = {
            "fastp": fastp_result.get("success", False),
            "star": star_result.get("success", False),
            "hisat2": hisat2_result.get("success", False),
            "featurecounts": fc_result.get("success", False)
        }

        # 收集解析错误
        parsing_errors = {}
        for step, result in [("fastp", fastp_result), ("star", star_result),
                           ("hisat2", hisat2_result), ("featurecounts", fc_result)]:
            if not result.get("success", False):
                parsing_errors[step] = result.get("error", "未知错误")

        # 样本数据对齐整合
        aligned_samples = _align_sample_metrics_integrated(
            fastp_result, star_result, hisat2_result, fc_result
        )

        return {
            "success": True,
            "results_directory": results_directory,
            "parsing_status": parsing_status,
            "parsing_errors": parsing_errors,
            "aligned_samples": aligned_samples,
            "summary": f"成功解析{sum(parsing_status.values())}/4个步骤，对齐{aligned_samples.get('total_count', 0)}个样本"
        }

    except Exception as e:
        logger.error(f"解析流水线结果失败: {e}")
        return {
            "success": False,
            "error": f"解析流水线结果失败: {str(e)}"
        }

def _align_sample_metrics_integrated(fastp_data, star_data, hisat2_data, fc_data):
    """整合版本的样本指标对齐函数"""
    # 从数据中收集所有样本ID
    all_sample_ids = set()

    for data_source in [fastp_data, star_data, hisat2_data, fc_data]:
        if data_source and data_source.get("success") and "sample_metrics" in data_source:
            for sample in data_source["sample_metrics"]:
                sid = sample.get("sample_id")
                if sid:
                    all_sample_ids.add(sid)

    # 构建样本指标映射
    fastp_samples = {}
    star_samples = {}
    hisat2_samples = {}
    fc_samples = {}
    
    if fastp_data and fastp_data.get("success"):
        fastp_samples = {s.get("sample_id"): s for s in fastp_data.get("sample_metrics", []) if s.get("sample_id")}
    
    if star_data and star_data.get("success"):
        star_samples = {s.get("sample_id"): s for s in star_data.get("sample_metrics", []) if s.get("sample_id")}
    
    if hisat2_data and hisat2_data.get("success"):
        hisat2_samples = {s.get("sample_id"): s for s in hisat2_data.get("sample_metrics", []) if s.get("sample_id")}
    
    if fc_data and fc_data.get("success"):
        fc_samples = {s.get("sample_id"): s for s in fc_data.get("sample_metrics", []) if s.get("sample_id")}

    # 对齐样本数据
    aligned_samples = []
    for sample_id in sorted(list(all_sample_ids)):
        sample_data = {
            "sample_id": sample_id,
            "fastp": fastp_samples.get(sample_id, {"error": "数据缺失"}),
            "star": star_samples.get(sample_id, {"error": "数据缺失"}),
            "hisat2": hisat2_samples.get(sample_id, {"error": "数据缺失"}),
            "featurecounts": fc_samples.get(sample_id, {"error": "数据缺失"}),
            "notes": []
        }
        
        # 添加简单的质量标记
        notes = []
        
        # FastP质量检查
        fastp_info = sample_data["fastp"]
        if not fastp_info.get("error"):
            q30_rate = fastp_info.get("q30_rate", 0)
            retention_rate = fastp_info.get("reads_passed_rate", 0)
            if q30_rate > 85 and retention_rate > 80:
                notes.append("FastP质量优秀")
            elif q30_rate < 70 or retention_rate < 60:
                notes.append("FastP质量需关注")
        
        # 比对质量检查（STAR或HISAT2）
        star_info = sample_data["star"]
        hisat2_info = sample_data["hisat2"]
        
        if not star_info.get("error"):
            unique_rate = star_info.get("uniquely_mapped_percentage", 0)
            if unique_rate > 90:
                notes.append("STAR比对优秀")
            elif unique_rate < 70:
                notes.append("STAR比对率偏低")
        elif not hisat2_info.get("error"):
            # 如果有HISAT2数据，检查其质量
            align_rate = hisat2_info.get("overall_alignment_rate", 0)
            if align_rate > 90:
                notes.append("HISAT2比对优秀")
            elif align_rate < 70:
                notes.append("HISAT2比对率偏低")
        
        # FeatureCounts质量检查
        fc_info = sample_data["featurecounts"]
        if not fc_info.get("error"):
            assign_rate = fc_info.get("assignment_rate", 0)
            if assign_rate > 0.8:
                notes.append("基因分配优秀")
            elif assign_rate < 0.6:
                notes.append("基因分配率偏低")
        
        sample_data["notes"] = notes
        aligned_samples.append(sample_data)

    return {
        "samples": aligned_samples,
        "total_count": len(aligned_samples),
        "summary": f"成功对齐{len(aligned_samples)}个样本的多步骤数据"
    }

















