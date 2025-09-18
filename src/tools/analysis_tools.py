"""
RNA-seq智能分析助手 - 分析和报告工具

包含：
- parse_pipeline_results: 解析流水线结果数据
- assess_sample_quality: 评估样本质量和健康度
- write_analysis_report: 生成综合分析报告
- write_analysis_markdown: 生成分析报告
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
    """解析RNA-seq流水线的所有步骤结果

    Args:
        results_directory: 结果目录路径

    Returns:
        Dict包含所有步骤的解析结果和样本数据
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

        return {
            "success": True,
            "results_directory": results_directory,
            "fastp_data": fastp_result if fastp_result.get("success") else {},
            "star_data": star_result if star_result.get("success") else {},
            "hisat2_data": hisat2_result if hisat2_result.get("success") else {},
            "featurecounts_data": fc_result if fc_result.get("success") else {},
            "parsing_status": parsing_status,
            "parsing_errors": parsing_errors,
            "summary": f"成功解析{sum(parsing_status.values())}/4个步骤"
        }

    except Exception as e:
        logger.error(f"解析流水线结果失败: {e}")
        return {
            "success": False,
            "error": f"解析流水线结果失败: {str(e)}"
        }


@tool
def assess_sample_quality(
    pipeline_data: Dict[str, Any],
    sample_groups: str = "[]"
) -> Dict[str, Any]:
    """基于流水线数据评估样本质量和健康度

    Args:
        pipeline_data: 从parse_pipeline_results获得的流水线数据
        sample_groups: 样本分组信息JSON字符串

    Returns:
        Dict包含样本质量评估和对齐后的指标数据
    """
    try:
        if not pipeline_data.get("success"):
            return {
                "success": False,
                "error": "流水线数据解析失败，无法进行质量评估"
            }

        # 解析样本分组
        import json
        try:
            sample_groups_list = json.loads(sample_groups) if sample_groups else []
        except:
            sample_groups_list = []

        fastp_data = pipeline_data.get("fastp_data", {})
        star_data = pipeline_data.get("star_data", {})
        hisat2_data = pipeline_data.get("hisat2_data", {})
        fc_data = pipeline_data.get("featurecounts_data", {})

        # 样本ID对齐和指标合并
        aligned_samples = _align_sample_metrics(
            sample_groups_list, fastp_data, star_data, hisat2_data, fc_data
        )

        return {
            "success": True,
            "aligned_samples": aligned_samples,
            "summary": f"成功对齐{len(aligned_samples.get('samples', []))}个样本的指标数据",
            "available_steps": [step for step, success in pipeline_data.get("parsing_status", {}).items() if success]
        }

    except Exception as e:
        logger.error(f"样本质量评估失败: {e}")
        return {
            "success": False,
            "error": f"样本质量评估失败: {str(e)}"
        }


@tool
def write_analysis_report(
    quality_assessment: Dict[str, Any],
    pipeline_config: Dict[str, Any],
    output_directory: str,
    report_title: str = "RNA-seq Analysis Report"
) -> Dict[str, Any]:
    """生成完整的RNA-seq分析报告

    Args:
        quality_assessment: 样本质量评估结果
        pipeline_config: 流水线配置信息
        output_directory: 输出目录
        report_title: 报告标题

    Returns:
        Dict包含报告生成结果和文件路径
    """
    try:
        output_path = Path(output_directory)
        output_path.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # 构建完整报告数据
        report_data = {
            "pipeline": {
                "steps": ["fastp", "star", "hisat2", "featurecounts"],
                "species": _resolve_species_from_genome(pipeline_config.get("genome_version", "unknown")),
                "genome_version": pipeline_config.get("genome_version", "unknown")
            },
            "context": {
                "results_dir": output_directory,
                "timestamp": timestamp,
                "sample_count": len(quality_assessment.get("aligned_samples", {}).get("samples", []))
            },
            "assessment": quality_assessment
        }

        # 保存JSON报告
        json_file = output_path / f"analysis_report_{timestamp}.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, ensure_ascii=False, indent=2)

        # 生成Markdown报告
        md_file = output_path / f"analysis_summary_{timestamp}.md"
        markdown_content = _generate_markdown_report(report_data)
        with open(md_file, 'w', encoding='utf-8') as f:
            f.write(markdown_content)

        return {
            "success": True,
            "json_report": str(json_file),
            "markdown_report": str(md_file),
            "report_data": report_data,
            "message": f"成功生成分析报告: {report_title}"
        }

    except Exception as e:
        logger.error(f"生成分析报告失败: {e}")
        return {
            "success": False,
            "error": f"生成分析报告失败: {str(e)}"
        }


def _resolve_species_from_genome(genome_version: str) -> str:
    """从基因组版本解析物种信息"""
    genome_to_species = {
        "hg38": "Homo sapiens (Human)",
        "hg19": "Homo sapiens (Human)",
        "mm39": "Mus musculus (Mouse)",
        "mm10": "Mus musculus (Mouse)",
        "mm9": "Mus musculus (Mouse)",
        "dm6": "Drosophila melanogaster (Fruit fly)",
        "danRer11": "Danio rerio (Zebrafish)",
        "xenLae2": "Xenopus laevis (African clawed frog)"
    }
    return genome_to_species.get(genome_version, "Unknown species")


def _align_sample_metrics(sample_groups, fastp_data, star_data, hisat2_data, fc_data):
    """对齐样本指标数据"""
    # 从数据中收集所有样本ID
    all_sample_ids = set()

    for data_source in [fastp_data, star_data, hisat2_data, fc_data]:
        if data_source and "sample_metrics" in data_source:
            for sample in data_source["sample_metrics"]:
                sid = sample.get("sample_id")
                if sid:
                    all_sample_ids.add(sid)

    # 如果有配置的样本组，优先使用
    expected_sample_ids = []
    if sample_groups:
        for group in sample_groups:
            sid = group.get("sample_id")
            if sid:
                expected_sample_ids.append(sid)
    else:
        expected_sample_ids = sorted(list(all_sample_ids))

    # 构建样本指标映射
    fastp_samples = {s.get("sample_id"): s for s in fastp_data.get("sample_metrics", []) if s.get("sample_id")}
    star_samples = {s.get("sample_id"): s for s in star_data.get("sample_metrics", []) if s.get("sample_id")}
    hisat2_samples = {s.get("sample_id"): s for s in hisat2_data.get("sample_metrics", []) if s.get("sample_id")}
    fc_samples = {s.get("sample_id"): s for s in fc_data.get("sample_metrics", []) if s.get("sample_id")}

    # 对齐样本数据
    aligned_samples = []
    for sample_id in expected_sample_ids:
        sample_data = {
            "sample_id": sample_id,
            "fastp": fastp_samples.get(sample_id, {"error": "数据缺失"}),
            "star": star_samples.get(sample_id, {"error": "数据缺失"}),
            "hisat2": hisat2_samples.get(sample_id, {"error": "数据缺失"}),
            "featurecounts": fc_samples.get(sample_id, {"error": "数据缺失"}),
            "notes": []
        }
        aligned_samples.append(sample_data)

    return {
        "samples": aligned_samples,
        "total_count": len(aligned_samples)
    }


def _generate_markdown_report(report_data: Dict[str, Any]) -> str:
    """生成Markdown格式报告"""
    pipeline = report_data.get("pipeline", {})
    context = report_data.get("context", {})
    assessment = report_data.get("assessment", {})

    lines = [
        f"# RNA-seq Analysis Report",
        "",
        f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Species**: {pipeline.get('species', 'Unknown')}",
        f"**Genome**: {pipeline.get('genome_version', 'Unknown')}",
        f"**Samples**: {context.get('sample_count', 0)}",
        "",
        "## Pipeline Steps",
        ""
    ]

    steps = pipeline.get("steps", [])
    for step in steps:
        lines.append(f"- {step.upper()}")

    lines.extend([
        "",
        "## Sample Summary",
        ""
    ])

    samples = assessment.get("aligned_samples", {}).get("samples", [])
    if samples:
        lines.extend([
            "| Sample ID | FastP | STAR | HISAT2 | FeatureCounts |",
            "|-----------|--------|------|--------|---------------|"
        ])

        for sample in samples:
            sample_id = sample.get("sample_id", "Unknown")
            fastp_status = "✅" if not sample.get("fastp", {}).get("error") else "❌"
            star_status = "✅" if not sample.get("star", {}).get("error") else "❌"
            hisat2_status = "✅" if not sample.get("hisat2", {}).get("error") else "❌"
            fc_status = "✅" if not sample.get("featurecounts", {}).get("error") else "❌"

            lines.append(f"| {sample_id} | {fastp_status} | {star_status} | {hisat2_status} | {fc_status} |")

    lines.extend([
        "",
        "---",
        "*This report was generated by RNA-seq智能分析助手*"
    ])

    return "\n".join(lines)


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
        
        # 获取metrics数据
        metrics = report.get("metrics", {})

        # 样本健康度摘要
        summary = report.get("summary", {})
        samples_info = summary.get("samples", {})
        if samples_info:
            lines.extend([
                "### 样本健康度评估",
                "",
                f"- **总样本数**: {samples_info.get('total', 'N/A')}",
                f"- **健康样本**: {samples_info.get('pass', 'N/A')}",
                f"- **异常样本**: {samples_info.get('fail', 'N/A')}",
                f"- **警告样本**: {samples_info.get('warn', 'N/A')}",
                f"- **整体状态**: {summary.get('status', 'N/A')}",
                "",
            ])

        # FastP质控结果
        fastp_metrics = metrics.get("fastp", {})
        if fastp_metrics:
            overall_stats = fastp_metrics.get("overall", {})
            lines.extend([
                "## FastP 质量控制",
                "",
                f"- **总输入reads**: {overall_stats.get('total_reads_before', 'N/A'):,}",
                f"- **总输出reads**: {overall_stats.get('total_reads_after', 'N/A'):,}",
                f"- **数据保留率**: {overall_stats.get('read_retention_rate', 'N/A')}%",
                f"- **平均Q20质量**: {overall_stats.get('average_q20_rate', 'N/A')}%",
                f"- **平均Q30质量**: {overall_stats.get('average_q30_rate', 'N/A')}%",
                f"- **平均GC含量**: {overall_stats.get('average_gc_content', 'N/A')}",
                "",
            ])

        # STAR比对结果
        star_metrics = metrics.get("star", {})
        if star_metrics and star_metrics.get("overall"):
            overall_stats = star_metrics.get("overall", {})
            lines.extend([
                "## STAR 序列比对",
                "",
                f"- **总输入reads**: {overall_stats.get('total_input_reads', 'N/A'):,}",
                f"- **唯一比对reads**: {overall_stats.get('total_uniquely_mapped', 'N/A'):,}",
                f"- **唯一比对率**: {overall_stats.get('overall_uniquely_mapped_percentage', 'N/A')}%",
                f"- **多重比对率**: {overall_stats.get('overall_multimapped_percentage', 'N/A')}%",
                f"- **总比对率**: {overall_stats.get('overall_mapped_percentage', 'N/A')}%",
                "",
            ])

        # HISAT2比对结果
        hisat2_metrics = metrics.get("hisat2", {})
        if hisat2_metrics and hisat2_metrics.get("overall"):
            overall_stats = hisat2_metrics.get("overall", {})
            lines.extend([
                "## HISAT2 序列比对",
                "",
                f"- **总输入reads**: {overall_stats.get('total_input_reads', 'N/A'):,}",
                f"- **唯一比对reads**: {overall_stats.get('total_uniquely_mapped', 'N/A'):,}",
                f"- **唯一比对率**: {overall_stats.get('overall_uniquely_mapped_percentage', 'N/A')}%",
                f"- **多重比对率**: {overall_stats.get('overall_multimapped_percentage', 'N/A')}%",
                f"- **总比对率**: {overall_stats.get('overall_mapped_percentage', 'N/A')}%",
                "",
            ])

        # FeatureCounts定量结果
        fc_metrics = metrics.get("featurecounts", {})
        if fc_metrics:
            overall_stats = fc_metrics.get("overall", {})
            lines.extend([
                "## FeatureCounts 基因定量",
                "",
                f"- **总输入reads**: {overall_stats.get('total_reads', 'N/A'):,}",
                f"- **分配到基因的reads**: {overall_stats.get('total_assigned', 'N/A'):,}",
                f"- **基因分配率**: {overall_stats.get('overall_assignment_rate', 0):.2%}",
                f"- **无特征reads**: {overall_stats.get('total_unassigned_nofeatures', 'N/A'):,}",
                f"- **多重映射reads**: {overall_stats.get('total_unassigned_ambiguity', 'N/A'):,}",
                f"- **质量过滤reads**: {overall_stats.get('total_unassigned_mappingquality', 'N/A'):,}",
                "",
            ])
        
        # 详细样本指标
        per_sample_data = report.get("per_sample", [])
        if per_sample_data:
            lines.extend([
                "## 详细样本指标",
                "",
                "| 样本ID | FastP保留率 | FeatureCounts分配率 | 健康状态 | 备注 |",
                "|--------|-------------|---------------------|----------|------|",
            ])

            for sample in per_sample_data:
                sample_id = sample.get("sample_id", "N/A")
                health_status = sample.get("health", "N/A")

                # FastP数据
                fastp_data = sample.get("fastp", {})
                fastp_rate = "N/A"
                if fastp_data and not fastp_data.get("error"):
                    fastp_rate = f"{fastp_data.get('reads_passed_rate', 'N/A'):.1f}"

                # FeatureCounts数据
                fc_data = sample.get("featurecounts", {})
                fc_rate = "N/A"
                if fc_data and not fc_data.get("error"):
                    fc_rate = f"{fc_data.get('assignment_rate', 0):.1%}"

                # 备注信息
                notes = sample.get("notes", [])
                notes_str = ", ".join(notes) if notes else "-"

                lines.append(f"| {sample_id} | {fastp_rate}% | {fc_rate} | {health_status} | {notes_str} |")

            lines.append("")
        
        # LLM分析内容
        if append_llm:
            llm_data = report.get("llm", {})
            if llm_data:
                lines.extend([
                    "## 智能分析与建议",
                    "",
                ])

                # 全局摘要
                global_summary = llm_data.get("global_summary", "")
                if global_summary:
                    lines.extend([
                        "### 全局分析摘要",
                        "",
                        global_summary,
                        "",
                    ])

                # 关键发现
                key_findings = llm_data.get("key_findings", [])
                if key_findings:
                    lines.extend([
                        "### 关键发现",
                        "",
                    ])
                    for finding in key_findings:
                        lines.append(f"- {finding}")
                    lines.append("")

                # 问题标记
                per_sample_flags = llm_data.get("per_sample_flags", [])
                if per_sample_flags:
                    lines.extend([
                        "### 样本问题标记",
                        "",
                    ])
                    for flag in per_sample_flags:
                        sample_id = flag.get("sample_id", "N/A")
                        issues = flag.get("issues", [])
                        severity = flag.get("severity", "unknown")
                        if issues:
                            lines.append(f"- **{sample_id}** ({severity}): {', '.join(issues)}")
                    lines.append("")

                # 改进建议
                recommendations = llm_data.get("recommendations", [])
                if recommendations:
                    lines.extend([
                        "### 改进建议",
                        "",
                    ])
                    for rec in recommendations:
                        action = rec.get("action", "")
                        reason = rec.get("reason", "")
                        if action:
                            lines.append(f"- **{action}**")
                            if reason:
                                lines.append(f"  - 原因: {reason}")
                    lines.append("")

                # 风险提示
                risks = llm_data.get("risks", [])
                if risks:
                    lines.extend([
                        "### 风险提示",
                        "",
                    ])
                    for risk in risks:
                        lines.append(f"- ⚠️ {risk}")
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