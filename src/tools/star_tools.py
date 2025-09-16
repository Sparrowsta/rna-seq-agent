"""
RNA-seq智能分析助手 - STAR比对工具

包含：
- run_nextflow_star: 执行STAR比对流程
- parse_star_metrics: 解析STAR比对结果
"""

import json
import subprocess
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块  
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.star")


@tool
def run_nextflow_star(
    star_params: Dict[str, Any],
    fastp_results: Dict[str, Any],
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None
) -> Dict[str, Any]:
    """执行 STAR 比对（精简版）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - STAR 索引优先使用 genome_info.star_index_dir；否则由 genome_info.fasta_path 推导
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 star.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行STAR比对"}

        fastp_results_dir = fastp_results.get("results_dir")
        if not fastp_results_dir:
            return {"success": False, "error": "FastP结果缺少results_dir"}

        per_sample = fastp_results.get("per_sample_outputs") or []
        if not per_sample:
            return {"success": False, "error": "FastP结果缺少per_sample_outputs"}

        # 2) 运行根目录与工作目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = Path(fastp_results_dir)
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"star_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"STAR启动：samples={len(per_sample)} results={results_dir} work={work_dir}")

        # 3) 解析 STAR 索引目录 - 直接基于fasta路径目录
        star_index_dir = ""
        if isinstance(genome_info, dict):
            star_index_dir = genome_info.get("star_index_dir") or genome_info.get("index_dir") or ""
            if not star_index_dir:
                fasta_path = genome_info.get("fasta_path") or genome_info.get("fasta")
                if fasta_path:
                    # 直接基于fasta路径的父目录构建索引目录
                    star_index_dir = str(Path(fasta_path).parent / "star_index")

        if not star_index_dir:
            return {"success": False, "error": "缺少STAR索引目录（genome_info.star_index_dir 或 fasta_path 必须提供）"}

        star_index_path = Path(star_index_dir)
        if not star_index_path.exists():
            logger.warning(f"STAR索引不存在: {star_index_path}")
            return {"success": False, "error": f"STAR索引不存在: {star_index_path}"}

        # 4) 构造 sample_inputs（仅使用 FastP 返回的结构）
        sample_inputs: List[Dict[str, Any]] = []
        for i, fp in enumerate(per_sample):
            sid = fp.get("sample_id", f"sample_{i+1}")
            r1 = fp.get("trimmed_single") or fp.get("trimmed_r1")
            r2 = fp.get("trimmed_r2")
            if not r1:
                continue
            sample_inputs.append({
                "sample_id": sid,
                "is_paired": bool(r2),
                "read1": r1,
                **({"read2": r2} if r2 else {})
            })
        if not sample_inputs:
            return {"success": False, "error": "未从FastP结果构造到任何样本输入"}

        # 5) 组装 Nextflow 参数
        cleaned_params: Dict[str, Any] = {}
        for k, v in (star_params or {}).items():
            if v is None or k in {"star_cpus", "outFileNamePrefix"}:
                continue
            cleaned_params[k.lstrip('-')] = v

        nf_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            "star_index": str(star_index_path),
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **cleaned_params,
        }

        # 保存参数文件到star子目录 - 改为版本化参数文件
        # M3: 参数版本化 - 将实际用于运行的nf_params写入版本化文件
        versioned_params_file = None
        try:
            # 构建临时 state 用于参数版本化
            from ..state import AgentState
            from .utils_tools import write_params_file
            temp_state = AgentState()
            temp_state.results_dir = str(results_dir)
            temp_state.results_timestamp = results_timestamp or ""
            
            # 写入版本化文件 - 使用实际运行的nf_params
            versioned_params_file = write_params_file("star", nf_params, temp_state)
            logger.info(f"STAR实际运行参数版本化文件已写入: {versioned_params_file}")
            
        except Exception as e:
            logger.warning(f"STAR参数版本化写入失败 (不影响执行): {e}")

        # 使用版本化文件作为参数文件，如果失败则降级到临时文件
        if versioned_params_file and versioned_params_file.exists():
            params_file = versioned_params_file
            logger.info(f"使用版本化参数文件运行STAR: {params_file}")
        else:
            # 降级方案：生成临时参数文件
            logger.warning("版本化参数文件不可用，使用临时参数文件")
            star_dir = results_dir / "star"
            star_dir.mkdir(parents=True, exist_ok=True)
            params_file = star_dir / "star_params.json"
            with open(params_file, "w", encoding="utf-8") as f:
                json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "star.nf",
            Path("/src/nextflow/star.nf"),
        ]
        nextflow_script = next((p for p in nf_candidates if p.exists()), None)
        if nextflow_script is None:
            return {"success": False, "error": "未找到 star.nf 脚本，请检查 /src/nextflow/star.nf 路径", "searched": [str(p) for p in nf_candidates]}

        logger.info(f"执行STAR比对 - 参数文件: {params_file}")
        logger.info(f"STAR索引: {nf_params['star_index']}")
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 star.nf publishDir 对齐）
        star_out = results_dir / "star"
        per_sample_outputs: List[Dict[str, Any]] = []
        for item in sample_inputs:
            sid = item["sample_id"]
            sdir = star_out / sid
            entry = {
                "sample_id": sid,
                "aligned_bam": str(sdir / f"{sid}.Aligned.sortedByCoord.out.bam"),
                "log_final": str(sdir / f"{sid}.Log.final.out"),
                "log_out": str(sdir / f"{sid}.Log.out"),
                "log_progress": str(sdir / f"{sid}.Log.progress.out"),
                "splice_junctions": str(sdir / f"{sid}.SJ.out.tab"),
            }
            qm = str(nf_params.get("quantMode", ""))
            if "TranscriptomeSAM" in qm:
                entry["transcriptome_bam"] = str(sdir / f"{sid}.Aligned.toTranscriptome.out.bam")
            if "GeneCounts" in qm:
                entry["gene_counts"] = str(sdir / f"{sid}.ReadsPerGene.out.tab")
            per_sample_outputs.append(entry)

        payload = {
            "success": result.returncode == 0,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": len(sample_inputs),
            "per_sample_outputs": per_sample_outputs,
        }
        if get_tools_config().settings.debug_mode:
            payload.update({"stdout": result.stdout, "stderr": result.stderr, "cmd": " ".join(cmd)})
        # 记录完成/失败日志
        if payload["success"]:
            logger.info(f"STAR完成：samples={len(sample_inputs)} results={results_dir}")
        else:
            logger.warning(f"STAR失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        return payload

    except Exception as e:
        logger.error(f"STAR异常：{e}")
        return {"success": False, "error": f"执行STAR比对失败: {str(e)}"}


@tool
def parse_star_metrics(results_directory: str) -> Dict[str, Any]:
    """解析STAR比对结果文件，提取比对指标
    
    Args:
        results_directory: STAR结果目录路径
    
    Returns:
        Dict: 解析后的比对指标数据
        {
            "success": bool,
            "total_samples": int,
            "results_directory": str,
            "per_sample_metrics": List[Dict],    # 每个样本的详细指标
            "overall_statistics": Dict,          # 总体统计信息
            "quality_assessment": Dict           # 质量评估
        }
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {
                "success": False,
                "error": f"STAR结果目录不存在: {results_directory}"
            }
        
        sample_metrics = []
        overall_stats = {
            "total_input_reads": 0,
            "total_mapped_reads": 0,
            "total_uniquely_mapped": 0,
            "total_multi_mapped": 0,
            "mapping_rates": [],
            "unique_mapping_rates": [],
            "multi_mapping_rates": [],
            "mismatch_rates": []
        }
        
        # 查找STAR输出目录
        star_dir = results_path / "star"
        if not star_dir.exists():
            try:
                logger.warning(f"STAR输出目录不存在: {star_dir}")
            except Exception:
                pass
            return {
                "success": False,
                "error": f"STAR输出目录不存在: {star_dir}"
            }
        
        # 遍历样本目录，查找Log.final.out文件
        for sample_dir in star_dir.iterdir():
            if not sample_dir.is_dir():
                continue
                
            # 与 star.nf 中 outFileNamePrefix 保持一致：{sample_id}/{sample_id}.*
            log_final_file = sample_dir / f"{sample_dir.name}.Log.final.out"
            if not log_final_file.exists():
                sample_metrics.append({
                    "sample_id": sample_dir.name,
                    "error": "Log.final.out文件不存在"
                })
                continue
            
            try:
                # 解析Log.final.out文件
                with open(log_final_file, 'r', encoding='utf-8') as f:
                    log_content = f.read()
                
                # 提取关键指标（使用正则表达式）
                input_reads = _extract_metric(log_content, r"Number of input reads.*?(\d+)")
                uniquely_mapped = _extract_metric(log_content, r"Uniquely mapped reads number.*?(\d+)")
                multi_mapped = _extract_metric(log_content, r"Number of reads mapped to multiple loci.*?(\d+)")
                unmapped = _extract_metric(log_content, r"Number of reads unmapped.*?(\d+)")
                
                # 计算比率
                if input_reads > 0:
                    mapping_rate = (uniquely_mapped + multi_mapped) / input_reads
                    unique_mapping_rate = uniquely_mapped / input_reads
                    multi_mapping_rate = multi_mapped / input_reads
                else:
                    mapping_rate = unique_mapping_rate = multi_mapping_rate = 0
                
                # 提取mismatch率
                mismatch_rate = _extract_metric(log_content, r"Mismatch rate per base.*?([\d.]+)%") / 100
                
                sample_metric = {
                    "sample_id": sample_dir.name,
                    "input_reads": input_reads,
                    "uniquely_mapped": uniquely_mapped,
                    "multi_mapped": multi_mapped,
                    "unmapped": unmapped,
                    "mapping_rate": round(mapping_rate, 4),
                    "unique_mapping_rate": round(unique_mapping_rate, 4),
                    "multi_mapping_rate": round(multi_mapping_rate, 4),
                    "mismatch_rate": round(mismatch_rate, 4),
                    "log_file": str(log_final_file)
                }
                
                sample_metrics.append(sample_metric)
                
                # 累加到总体统计
                overall_stats["total_input_reads"] += input_reads
                overall_stats["total_mapped_reads"] += (uniquely_mapped + multi_mapped)
                overall_stats["total_uniquely_mapped"] += uniquely_mapped
                overall_stats["total_multi_mapped"] += multi_mapped
                overall_stats["mapping_rates"].append(mapping_rate)
                overall_stats["unique_mapping_rates"].append(unique_mapping_rate)
                overall_stats["multi_mapping_rates"].append(multi_mapping_rate)
                overall_stats["mismatch_rates"].append(mismatch_rate)
                
            except Exception as e:
                sample_metrics.append({
                    "sample_id": sample_dir.name,
                    "log_file": str(log_final_file),
                    "error": f"解析失败: {str(e)}"
                })
        
        # 计算总体指标
        total_samples = len([m for m in sample_metrics if "error" not in m])
        if total_samples > 0:
            overall_mapping_rate = overall_stats["total_mapped_reads"] / overall_stats["total_input_reads"]
            overall_unique_rate = overall_stats["total_uniquely_mapped"] / overall_stats["total_input_reads"]
            overall_multi_rate = overall_stats["total_multi_mapped"] / overall_stats["total_input_reads"]
            avg_mismatch_rate = sum(overall_stats["mismatch_rates"]) / len(overall_stats["mismatch_rates"])
        else:
            overall_mapping_rate = overall_unique_rate = overall_multi_rate = avg_mismatch_rate = 0
        
        # 质量评估
        quality_assessment = {
            "overall_quality": "good" if overall_mapping_rate > 0.85 else "moderate" if overall_mapping_rate > 0.7 else "poor",
            "unique_mapping_status": "good" if overall_unique_rate > 0.8 else "moderate" if overall_unique_rate > 0.6 else "poor",
            "multi_mapping_status": "good" if overall_multi_rate < 0.2 else "moderate" if overall_multi_rate < 0.3 else "high"
        }
        
        result = {
            "success": True,
            "total_samples": total_samples,
            "results_directory": results_directory,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_input_reads": overall_stats["total_input_reads"],
                "total_mapped_reads": overall_stats["total_mapped_reads"],
                "total_uniquely_mapped": overall_stats["total_uniquely_mapped"],
                "total_multi_mapped": overall_stats["total_multi_mapped"],
                "overall_mapping_rate": round(overall_mapping_rate, 4),
                "overall_unique_mapping_rate": round(overall_unique_rate, 4),
                "overall_multi_mapping_rate": round(overall_multi_rate, 4),
                "average_mismatch_rate": round(avg_mismatch_rate, 4)
            },
            "quality_assessment": quality_assessment
        }
        try:
            logger.info(
                f"STAR结果: samples={result['total_samples']} map={result['overall_statistics']['overall_mapping_rate']} "
                f"unique={result['overall_statistics']['overall_unique_mapping_rate']} multi={result['overall_statistics']['overall_multi_mapping_rate']}"
            )
            if sample_metrics:
                logger.debug(f"STAR样本预览：{sample_metrics[0]}")
        except Exception:
            pass
        return result
        
    except Exception as e:
        try:
            logger.error(f"解析STAR结果失败：{e}")
        except Exception:
            pass
        return {
            "success": False,
            "error": f"解析STAR结果失败: {str(e)}"
        }


def _extract_metric(text: str, pattern: str) -> float:
    """从文本中提取数值指标的辅助函数"""
    match = re.search(pattern, text)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return 0.0
    return 0.0