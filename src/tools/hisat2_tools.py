"""
RNA-seq智能分析助手 - HISAT2工具集

包含：
- run_nextflow_hisat2: 执行HISAT2比对
- parse_hisat2_metrics: 解析HISAT2结果指标
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
from ..config import get_tools_config, Settings
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.hisat2")


@tool
def run_nextflow_hisat2(
    hisat2_params: Dict[str, Any],
    fastp_results: Dict[str, Any],
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None
) -> Dict[str, Any]:
    """执行 HISAT2 比对（精简版，与 run_nextflow_star 等价）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - HISAT2 索引优先使用 genome_info.hisat2_index_dir；否则由 genome_info.fasta_path 推导
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 hisat2.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # M3: 参数版本化 - 在实际执行前写入参数版本化文件
        try:
            # 从 fastp_results 获取 results_dir，直接传递给write_params_file
            fastp_results_dir = fastp_results.get("results_dir") if fastp_results else ""
            if fastp_results_dir:
                from .utils_tools import write_params_file

                # 直接传递results_dir，不需要创建临时state
                param_file_path = write_params_file(
                    "hisat2",
                    hisat2_params,
                    results_dir=fastp_results_dir
                )
                logger.info(f"HISAT2执行参数版本化文件已写入: {param_file_path}")

        except Exception as e:
            logger.warning(f"HISAT2参数版本化写入失败 (不影响执行): {e}")

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行HISAT2比对"}

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
        work_dir = tools_config.settings.data_dir / "work" / f"hisat2_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)

        # 3) 解析 HISAT2 索引目录 - 直接基于fasta路径目录
        hisat2_index_prefix = ""
        if isinstance(genome_info, dict):
            hisat2_index_dir = genome_info.get("hisat2_index_dir") or genome_info.get("index_dir") or ""
            if hisat2_index_dir:
                hisat2_index_prefix = str(Path(hisat2_index_dir) / "genome")
            else:
                fasta_path_raw = genome_info.get("fasta_path") or genome_info.get("fasta")
                if fasta_path_raw:
                    fasta_path = fasta_path_raw
                    # 直接基于fasta路径的父目录构建索引目录
                    hisat2_index_prefix = str(Path(fasta_path).parent / "hisat2_index" / "genome")

        if not hisat2_index_prefix:
            return {"success": False, "error": "缺少HISAT2索引目录（genome_info.hisat2_index_dir 或 fasta_path 必须提供）"}

        # 检查索引文件是否存在
        index_files = list(Path(hisat2_index_prefix).parent.glob(f"{Path(hisat2_index_prefix).name}.*.ht2"))
        if not index_files:
            return {"success": False, "error": f"HISAT2索引文件不存在: {hisat2_index_prefix}.*.ht2"}

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
        for k, v in (hisat2_params or {}).items():
            if v is None or k in {"hisat2_cpus", "threads", "p"}:
                continue
            cleaned_params[k.lstrip('-')] = v

        nf_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            "hisat2_index": hisat2_index_prefix,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **cleaned_params,
        }

        # 参数版本化 - 使用新的接口直接传递results_dir
        versioned_params_file = None
        try:
            from .utils_tools import write_params_file
            # 直接传递results_dir，不需要创建临时state
            versioned_params_file = write_params_file(
                "hisat2",
                nf_params,
                results_dir=str(results_dir)
            )
            logger.info(f"HISAT2实际运行参数版本化文件已写入: {versioned_params_file}")

        except Exception as e:
            logger.warning(f"HISAT2参数版本化写入失败 (不影响执行): {e}")

        # 使用版本化文件作为参数文件，如果失败则降级到临时文件
        if versioned_params_file and versioned_params_file.exists():
            params_file = versioned_params_file
            logger.info(f"使用版本化参数文件运行HISAT2: {params_file}")
        else:
            # 降级方案：生成临时参数文件
            logger.warning("版本化参数文件不可用，使用临时参数文件")
            hisat2_dir = results_dir / "hisat2"
            hisat2_dir.mkdir(parents=True, exist_ok=True)
            params_file = hisat2_dir / "hisat2_params.json"
            with open(params_file, "w", encoding="utf-8") as f:
                json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "hisat2.nf"
        if not nextflow_script.exists():
            return {
                "success": False, 
                "error": "未找到 hisat2.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "hisat2.nf")]
            }

        logger.info(f"执行HISAT2比对 - 参数文件: {params_file}")
        logger.info(f"HISAT2索引: {nf_params['hisat2_index']}")
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 hisat2.nf publishDir 对齐）
        hisat2_out = results_dir / "hisat2"
        per_sample_outputs: List[Dict[str, Any]] = []
        for item in sample_inputs:
            sid = item["sample_id"]
            sdir = hisat2_out / sid
            entry = {
                "sample_id": sid,
                "aligned_bam": str(sdir / f"{sid}.hisat2.bam"),
                "align_summary": str(sdir / f"{sid}.align_summary.txt"),
                "bam_index": str(sdir / f"{sid}.hisat2.bam.bai"),
            }
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
            payload.update({"stderr": result.stderr, "cmd": " ".join(cmd)})
        return payload

    except Exception as e:
        return {"success": False, "error": f"执行HISAT2比对失败: {str(e)}"}


@tool
def parse_hisat2_metrics(results_directory: str) -> Dict[str, Any]:
    """解析HISAT2比对结果文件，提取比对指标
    
    Args:
        results_directory: HISAT2结果目录路径
        
    Returns:
        Dict: 包含样本指标和整体统计的结果
    """
    try:
        results_path = Path(results_directory)
        
        if not results_path.exists():
            return {"success": False, "error": f"结果目录不存在: {results_directory}"}
        
        # 查找HISAT2结果子目录
        hisat2_dir = results_path / "hisat2"
        if not hisat2_dir.exists():
            return {"success": False, "error": f"HISAT2子目录不存在: {hisat2_dir}"}
        
        sample_metrics = []
        overall_stats = {
            "total_input_reads": 0,
            "total_mapped_reads": 0,
            "total_uniquely_mapped": 0,
            "total_multi_mapped": 0,
            "mapping_rates": [],
            "unique_mapping_rates": [],
            "multi_mapping_rates": []
        }
        
        # 扫描样本目录
        sample_dirs = [d for d in hisat2_dir.iterdir() if d.is_dir()]
        total_samples = len(sample_dirs)
        
        if total_samples == 0:
            return {"success": False, "error": "未找到任何样本目录"}
        
        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            summary_file = sample_dir / f"{sample_id}.align_summary.txt"
            
            if not summary_file.exists():
                continue
            
            # 解析HISAT2比对摘要
            with open(summary_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # HISAT2输出格式解析
            metrics = _parse_hisat2_summary(content, sample_id)
            if metrics:
                sample_metrics.append(metrics)
                
                # 累积整体统计
                overall_stats["total_input_reads"] += metrics["input_reads"]
                overall_stats["total_mapped_reads"] += metrics["mapped_reads"]
                overall_stats["total_uniquely_mapped"] += metrics.get("uniquely_mapped", 0)
                overall_stats["total_multi_mapped"] += metrics.get("multi_mapped", 0)
                overall_stats["mapping_rates"].append(metrics["mapping_rate"])
                overall_stats["unique_mapping_rates"].append(metrics.get("unique_mapping_rate", 0))
                overall_stats["multi_mapping_rates"].append(metrics.get("multi_mapping_rate", 0))
        
        # 计算整体比率
        if overall_stats["total_input_reads"] > 0:
            overall_mapping_rate = overall_stats["total_mapped_reads"] / overall_stats["total_input_reads"]
            overall_unique_rate = overall_stats["total_uniquely_mapped"] / overall_stats["total_input_reads"]
            overall_multi_rate = overall_stats["total_multi_mapped"] / overall_stats["total_input_reads"]
        else:
            overall_mapping_rate = overall_unique_rate = overall_multi_rate = 0
        
        # 质量评估
        quality_assessment = {
            "overall_quality": "good" if overall_mapping_rate > 0.85 else "moderate" if overall_mapping_rate > 0.7 else "poor",
            "unique_mapping_status": "good" if overall_unique_rate > 0.8 else "moderate" if overall_unique_rate > 0.6 else "poor",
            "multi_mapping_status": "good" if overall_multi_rate < 0.2 else "moderate" if overall_multi_rate < 0.3 else "high"
        }
        
        return {
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
            },
            "quality_assessment": quality_assessment
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"解析HISAT2结果失败: {str(e)}"
        }


def _parse_hisat2_summary(content: str, sample_id: str) -> Dict[str, Any]:
    """解析HISAT2比对摘要内容（健壮版）

    兼容数字包含千分位逗号、单端/双端输出差异，支持从 overall alignment rate 回退计算。
    """
    try:
        lines = [ln.strip() for ln in (content or "").splitlines() if ln.strip()]
        metrics: Dict[str, Any] = {"sample_id": sample_id}

        # 工具函数：数字提取（支持1,234,567）
        def _parse_int_token(tok: str) -> int:
            try:
                return int(str(tok).replace(",", "").strip())
            except Exception:
                return 0

        # 提取总reads
        # 例："25,000,000 reads; of these:" 或 "25000000 reads; of these:"
        for ln in lines:
            if "reads; of these:" in ln:
                # 优先用正则
                m = re.search(r"([\d,]+)\s+reads; of these:", ln)
                if m:
                    metrics["input_reads"] = _parse_int_token(m.group(1))
                else:
                    # 回退：取首个token
                    parts = ln.split()
                    if parts:
                        metrics["input_reads"] = _parse_int_token(parts[0])
                break

        # 提取映射计数（双端优先匹配concordantly，单端匹配非concordantly）
        for ln in lines:
            if "aligned concordantly exactly 1 time" in ln:
                metrics["uniquely_mapped"] = _extract_reads_count(ln)
            elif "aligned concordantly >1 times" in ln:
                metrics["multi_mapped"] = _extract_reads_count(ln)
            elif ("aligned exactly 1 time" in ln) and ("concordantly" not in ln):
                metrics["uniquely_mapped"] = _extract_reads_count(ln)
            elif ("aligned >1 times" in ln) and ("concordantly" not in ln):
                metrics["multi_mapped"] = _extract_reads_count(ln)

        # 解析 overall alignment rate（百分比），作为映射率的兜底
        overall_rate = None
        for ln in lines[::-1]:
            # 例："90.00% overall alignment rate"
            if "overall alignment rate" in ln:
                m = re.search(r"([\d.]+)%\s+overall alignment rate", ln)
                if m:
                    try:
                        overall_rate = float(m.group(1)) / 100.0
                    except Exception:
                        overall_rate = None
                break

        # 计算派生指标
        if "input_reads" in metrics:
            total_reads = int(metrics.get("input_reads", 0) or 0)
            uniquely_mapped = int(metrics.get("uniquely_mapped", 0) or 0)
            multi_mapped = int(metrics.get("multi_mapped", 0) or 0)
            mapped_reads = uniquely_mapped + multi_mapped

            # 若unique+multi均未捕获，但有overall_rate，回退计算mapped_reads
            if mapped_reads == 0 and overall_rate is not None and total_reads > 0:
                mapped_reads = int(round(total_reads * overall_rate))

            unmapped_reads = max(0, total_reads - mapped_reads)
            metrics["mapped_reads"] = mapped_reads
            metrics["unmapped_reads"] = unmapped_reads

            if total_reads > 0:
                metrics["mapping_rate"] = round(mapped_reads / total_reads, 4)
                metrics["unique_mapping_rate"] = round(uniquely_mapped / total_reads, 4) if uniquely_mapped else 0.0
                metrics["multi_mapping_rate"] = round(multi_mapped / total_reads, 4) if multi_mapped else 0.0
                metrics["unmapped_rate"] = round(unmapped_reads / total_reads, 4)
            else:
                metrics.update({
                    "mapping_rate": 0.0,
                    "unique_mapping_rate": 0.0,
                    "multi_mapping_rate": 0.0,
                    "unmapped_rate": 0.0,
                })

        return metrics if "input_reads" in metrics else None

    except Exception as e:
        logger.error(f"解析HISAT2摘要失败 (样本 {sample_id}): {e}")
        return None


def _extract_reads_count(line: str) -> int:
    """从HISAT2输出行中提取reads数量（支持千分位逗号）"""
    match = re.search(r'([\d,]+)\s+\([\d.]+%\)', line.strip())
    if match:
        try:
            return int(match.group(1).replace(',', ''))
        except Exception:
            return 0
    return 0