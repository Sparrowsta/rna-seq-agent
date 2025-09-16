import json
from typing import Dict, Any
from datetime import datetime
from pathlib import Path
from ..state import AgentState
from ..tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files,
)
from ..config.settings import Settings
from ..logging_bootstrap import get_logger

logger = get_logger("rna.nodes.detect")


async def detect_node(state: AgentState) -> Dict[str, Any]:
    """Detect节点 - 直接调用工具执行全面检测（不依赖Plan）"""
    logger.info("正在执行全面环境与数据检测")

    # 生成时间戳和结果目录
    settings = Settings()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_data_path = str(settings.data_dir)
    results_dir = f"{base_data_path}/results/{timestamp}"
    
    logger.info(f"生成结果目录: {results_dir}")

    results: Dict[str, Any] = {}
    errors = []

    try:
        # 使用 BaseTool.invoke 调用以避免弃用警告
        results["analyze_fastq_data"] = scan_fastq_files.invoke({})
    except Exception as e:
        errors.append(f"scan_fastq_files: {e}")

    try:
        results["assess_system_readiness"] = scan_system_resources.invoke({})
    except Exception as e:
        errors.append(f"scan_system_resources: {e}")

    try:
        # 若用户之前指定过 genome_id，可传入聚焦
        user_req = getattr(state, "user_requirements", {}) or {}
        focus = None
        raw = (user_req.get("raw_input") or "").strip()
        for gid in ["hg38", "hg19", "mm39", "mm10", "rn6", "ce11", "dm6"]:
            if gid in raw.lower():
                focus = gid
                break

        # 首次调用scan_genome_files
        if focus:
            genome_result = scan_genome_files.invoke({"genome_id": focus})
        else:
            genome_result = scan_genome_files.invoke({})

        # 如果找不到genomes.json文件，创建一个空的
        if genome_result.get("detection_status") == "no_config":
            missing_config_file = genome_result.get("missing_config_file")
            if missing_config_file:
                try:
                    # 创建目录（如果不存在）
                    config_path = Path(missing_config_file)
                    config_path.parent.mkdir(parents=True, exist_ok=True)

                    # 创建空的genomes.json文件
                    with open(config_path, 'w', encoding='utf-8') as f:
                        json.dump({}, f, indent=2)

                    logger.info(f"已创建空的基因组配置文件: {config_path}")

                    # 设置创建后的结果状态 (适配新的scan_genome_files结构)
                    genome_result = {
                        "detection_status": "success",
                        "available_count": [],
                        "total_configured": 0,
                        "available_star_index": [],
                        "available_hisat2_index": []
                    }
                except Exception as create_error:
                    logger.warning(f"创建基因组配置文件失败: {create_error}")

        results["verify_genome_setup"] = genome_result
    except Exception as e:
        errors.append(f"scan_genome_files: {e}")

    # 工具可用性 - Docker环境保证所有工具可用
    for tool in ["fastp", "star", "hisat2", "featurecounts"]:
        results[f"check_{tool}_availability"] = {
            "tool_name": tool,
            "available": True,
            "environment": "docker",
            "guaranteed_by_docker": True
        }

    # 读取FASTQ统计，增强可观测性
    analyze_fastq_data = results.get("analyze_fastq_data") or {}
    fastq_total_samples = analyze_fastq_data.get("total_samples")
    if fastq_total_samples is None:
        fastq_total_samples = len((analyze_fastq_data.get("samples") or {}))
    fastq_total_files = analyze_fastq_data.get("total_files") or 0
    search_roots = ",".join(analyze_fastq_data.get("search_roots") or [])

    if search_roots:
        sample_names = list((analyze_fastq_data.get("samples") or {}).keys())
        preview = ", ".join(sample_names[:3]) + ("..." if len(sample_names) > 3 else "")
        logger.info(f"FASTQ扫描: roots=[{search_roots}] files={fastq_total_files} samples={fastq_total_samples} preview=[{preview}]")

    # 汇总
    # 获取可用基因组信息 (适配新的scan_genome_files结构)
    genome_setup = results.get('verify_genome_setup') or {}

    # 计算有索引的基因组总数 (STAR或HISAT2任一可用)
    star_genomes = set(genome_setup.get('available_star_index', []))
    hisat2_genomes = set(genome_setup.get('available_hisat2_index', []))
    available_indexed_genomes = star_genomes | hisat2_genomes

    summary_parts = [
        "检测完成",
        f"FASTQ样本: {fastq_total_samples}",
        f"可用基因组: {len(available_indexed_genomes)}",
        "工具: " + ", ".join(
            f"{t}:✅" for t in ["fastp", "star", "hisat2", "featurecounts"]
        ) + " (Docker保证)",
    ]
    if errors:
        summary_parts.append(f"错误 {len(errors)}")

    query_summary = " | ".join(summary_parts)
    logger.info(f"检测完成: {query_summary}")

    return {
        "success": True,
        "query_summary": query_summary,
        "status": "prepare",
        "query_results": results,
        "execution_errors": errors or None,
        "response": query_summary,
        
        # 时间戳和目录信息
        "results_dir": results_dir,
        "results_timestamp": timestamp,
        "base_data_path": base_data_path,
        
        # 为后续节点预设 nextflow_config
        "nextflow_config": {
            "results_dir": results_dir
        }
    }
