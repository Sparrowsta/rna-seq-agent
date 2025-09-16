"""
RNA-seq智能分析助手 - 分析和报告工具

包含：
- write_analysis_markdown: 生成分析报告
"""

import json
import time
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
def write_analysis_markdown(
    analysis_report: Dict[str, Any],
    output_path: str = ""
) -> Dict[str, Any]:
    """生成RNA-seq分析的Markdown报告
    
    Args:
        analysis_report: 分析报告数据，包含FastP、STAR、FeatureCounts结果
        output_path: 输出路径，默认自动生成
    
    Returns:
        生成结果字典，包含文件路径等信息
    """
    try:
        if not analysis_report or not isinstance(analysis_report, dict):
            return {
                "success": False,
                "error": "分析报告数据为空或格式错误"
            }
        
        # 确定输出路径
        if output_path:
            output_file = Path(output_path)
        else:
            # 从report中获取目录信息
            context = analysis_report.get("context", {})
            results_dir = context.get("results_dir")
            timestamp = context.get("timestamp")
            
            if results_dir and timestamp:
                target_dir = Path(results_dir) / "reports" / timestamp
                target_dir.mkdir(parents=True, exist_ok=True)
                output_file = target_dir / "analysis_summary.md"
            else:
                # 默认路径
                config = get_tools_config()
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                default_dir = config.settings.data_dir / "reports" / timestamp
                default_dir.mkdir(parents=True, exist_ok=True)
                output_file = default_dir / "analysis_summary.md"
        
        # 生成Markdown内容
        markdown_content = _render_analysis_markdown(analysis_report)
        
        # 写入文件
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(markdown_content)
        
        logger.info(f"分析报告生成完成: {output_file}")
        
        return {
            "success": True,
            "output_file": str(output_file),
            "file_size": output_file.stat().st_size,
            "message": f"分析报告已生成: {output_file.name}"
        }
        
    except Exception as e:
        logger.error(f"生成分析报告失败: {e}")
        return {
            "success": False,
            "error": f"生成分析报告失败: {str(e)}"
        }


def _render_analysis_markdown(report: Dict[str, Any], append_llm: bool = True) -> str:
    """渲染分析报告为Markdown格式
    
    Args:
        report: 分析报告数据
        append_llm: 是否添加LLM分析内容
    
    Returns:
        Markdown格式的报告内容
    """
    try:
        lines = []
        
        # 标题和基本信息
        context = report.get("context", {})
        timestamp = context.get("timestamp") or datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        lines.extend([
            "# RNA-seq 分析报告",
            "",
            f"**生成时间**: {timestamp}",
            "",
            "## 概览",
            "",
        ])
        
        # 样本健康度摘要
        health_summary = report.get("sample_health_summary", {})
        if health_summary:
            lines.extend([
                "### 样本健康度评估",
                "",
                f"- **总样本数**: {health_summary.get('total_samples', 'N/A')}",
                f"- **健康样本**: {health_summary.get('healthy_samples', 'N/A')}",
                f"- **异常样本**: {health_summary.get('problematic_samples', 'N/A')}",
                f"- **整体状态**: {health_summary.get('overall_status', 'N/A')}",
                "",
            ])
        
        # FastP质控结果
        fastp_summary = report.get("fastp_summary", {})
        if fastp_summary:
            overall_stats = fastp_summary.get("overall_statistics", {})
            lines.extend([
                "## FastP 质量控制",
                "",
                f"- **处理样本数**: {fastp_summary.get('sample_count', 'N/A')}",
                f"- **总输入reads**: {overall_stats.get('total_reads_before', 'N/A'):,}",
                f"- **总输出reads**: {overall_stats.get('total_reads_after', 'N/A'):,}",
                f"- **数据保留率**: {overall_stats.get('read_retention_rate', 'N/A')}%",
                f"- **平均Q30质量**: {overall_stats.get('average_q30_rate', 'N/A')}%",
                f"- **平均GC含量**: {overall_stats.get('average_gc_content', 'N/A')}%",
                "",
            ])
        
        # STAR比对结果
        star_summary = report.get("star_summary", {})
        if star_summary:
            overall_stats = star_summary.get("overall_statistics", {})
            lines.extend([
                "## STAR 序列比对",
                "",
                f"- **比对样本数**: {star_summary.get('sample_count', 'N/A')}",
                f"- **总输入reads**: {overall_stats.get('total_input_reads', 'N/A'):,}",
                f"- **唯一比对reads**: {overall_stats.get('total_uniquely_mapped', 'N/A'):,}",
                f"- **唯一比对率**: {overall_stats.get('overall_uniquely_mapped_percentage', 'N/A')}%",
                f"- **多重比对率**: {overall_stats.get('overall_multimapped_percentage', 'N/A')}%",
                f"- **总比对率**: {overall_stats.get('overall_mapped_percentage', 'N/A')}%",
                "",
            ])
        
        # FeatureCounts定量结果
        fc_summary = report.get("featurecounts_summary", {})
        if fc_summary:
            overall_stats = fc_summary.get("overall_statistics", {})
            lines.extend([
                "## FeatureCounts 基因定量",
                "",
                f"- **定量样本数**: {fc_summary.get('sample_count', 'N/A')}",
                f"- **总输入reads**: {overall_stats.get('total_reads', 'N/A'):,}",
                f"- **分配到基因的reads**: {overall_stats.get('assigned_reads', 'N/A'):,}",
                f"- **基因分配率**: {overall_stats.get('overall_assignment_rate', 'N/A')}%",
                f"- **无特征reads率**: {overall_stats.get('unassigned_no_features_rate', 'N/A')}%",
                f"- **多重映射reads率**: {overall_stats.get('unassigned_multimapping_rate', 'N/A')}%",
                "",
            ])
        
        # 详细样本指标
        sample_details = report.get("detailed_sample_metrics", [])
        if sample_details:
            lines.extend([
                "## 详细样本指标",
                "",
                "| 样本ID | FastP保留率 | STAR唯一比对率 | STAR总比对率 | FeatureCounts分配率 | 健康状态 |",
                "|--------|-------------|----------------|--------------|---------------------|----------|",
            ])
            
            for sample in sample_details:
                sample_id = sample.get("sample_id", "N/A")
                fastp_retention = sample.get("fastp_retention_rate", "N/A")
                star_unique = sample.get("star_unique_rate", "N/A")
                star_total = sample.get("star_total_rate", "N/A")
                fc_assignment = sample.get("fc_assignment_rate", "N/A")
                health_status = sample.get("health_status", "N/A")
                
                lines.append(f"| {sample_id} | {fastp_retention}% | {star_unique}% | {star_total}% | {fc_assignment}% | {health_status} |")
            
            lines.append("")
        
        # LLM分析内容
        if append_llm:
            llm_analysis = report.get("llm_analysis", {})
            if llm_analysis:
                lines.extend([
                    "## 智能分析与建议",
                    "",
                ])
                
                # 质量评估
                quality_assessment = llm_analysis.get("quality_assessment", "")
                if quality_assessment:
                    lines.extend([
                        "### 数据质量评估",
                        "",
                        quality_assessment,
                        "",
                    ])
                
                # 问题诊断
                issues_identified = llm_analysis.get("issues_identified", [])
                if issues_identified:
                    lines.extend([
                        "### 识别的问题",
                        "",
                    ])
                    for issue in issues_identified:
                        lines.append(f"- {issue}")
                    lines.append("")
                
                # 优化建议
                optimization_suggestions = llm_analysis.get("optimization_suggestions", [])
                if optimization_suggestions:
                    lines.extend([
                        "### 优化建议",
                        "",
                    ])
                    for suggestion in optimization_suggestions:
                        lines.append(f"- {suggestion}")
                    lines.append("")
        
        # 生成时间和版本信息
        lines.extend([
            "---",
            "",
            f"*报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
            "",
            "*由 RNA-seq智能分析助手 生成*",
        ])
        
        return "\n".join(lines)
        
    except Exception as e:
        logger.error(f"渲染Markdown报告失败: {e}")
        return f"# 报告生成失败\n\n错误: {str(e)}"