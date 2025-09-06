from typing import Dict, Any
from datetime import datetime
from ..state import AgentState
from ..tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files,
    check_tool_availability,
)
from ..config.settings import Settings


async def detect_node(state: AgentState) -> Dict[str, Any]:
    """Detect节点 - 直接调用工具执行全面检测（不依赖Plan）"""
    print("🔎 正在执行全面环境与数据检测…")

    # 生成时间戳和结果目录
    settings = Settings()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_data_path = str(settings.data_dir)
    results_dir = f"{base_data_path}/results/{timestamp}"
    
    print(f"📁 生成结果目录: {results_dir}")

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
        if focus:
            results["verify_genome_setup"] = scan_genome_files.invoke({"genome_id": focus})
        else:
            results["verify_genome_setup"] = scan_genome_files.invoke({})
    except Exception as e:
        errors.append(f"scan_genome_files: {e}")

    # 工具可用性
    for tool in ["fastp", "star", "hisat2", "featurecounts"]:
        try:
            results[f"check_{tool}_availability"] = check_tool_availability.invoke({"tool_name": tool})
        except Exception as e:
            errors.append(f"check_{tool}_availability: {e}")

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
        print(f"🔍 FASTQ扫描: roots=[{search_roots}] files={fastq_total_files} samples={fastq_total_samples} preview=[{preview}]")

    # 汇总
    summary_parts = [
        "检测完成",
        f"FASTQ样本: {fastq_total_samples}",
        f"可用基因组: {(results.get('verify_genome_setup') or {}).get('available_genomes', 0)}",
        "工具: " + ", ".join(
            f"{t}:{'✅' if (results.get(f'check_{t}_availability') or {}).get('available') else '❌'}"
            for t in ["fastp", "star", "hisat2", "featurecounts"]
        ),
    ]
    if errors:
        summary_parts.append(f"错误 {len(errors)}")

    query_summary = " | ".join(summary_parts)
    print(f"✅ 检测完成: {query_summary}")

    return {
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
