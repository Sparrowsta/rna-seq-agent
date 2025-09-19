"""
RNA-seq智能分析助手 - STAR比对工具

包含：
- run_nextflow_star: 执行STAR比对流程
- parse_star_metrics: 解析STAR比对结果
"""

import json
import subprocess
import re
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
    genome_paths: Dict[str, str],
    resource_config: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """执行 STAR 比对（精简版）

    Args:
        star_params: STAR执行参数
        fastp_results: FastP质控结果
        genome_paths: 简化的基因组路径信息（从节点动态提取）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - 直接使用 genome_paths["index_path"] 获取STAR索引路径
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 star.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 直接从传入的路径信息中提取所需字段
        genome_id = genome_paths.get("genome_id", "unknown")
        star_index_path = genome_paths.get("index_path", "")

        if not star_index_path:
            return {
                "success": False,
                "error": f"基因组路径信息缺少STAR索引路径: {genome_id}"
            }

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行STAR比对"}

        fastp_results_dir = fastp_results.get("results_dir")
        if not fastp_results_dir:
            return {"success": False, "error": "FastP结果缺少results_dir"}

        fastp_per_sample_outputs = fastp_results.get("per_sample_outputs") or []
        if not fastp_per_sample_outputs:
            return {"success": False, "error": "FastP结果缺少per_sample_outputs"}

        # 2) 运行根目录与工作目录
        results_dir = Path(fastp_results_dir)
        run_id = results_dir.name
        work_dir = tools_config.settings.data_dir / "work" / f"star_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"STAR启动:samples={len(fastp_per_sample_outputs)} results={results_dir} work={work_dir}")

        # 3) 直接使用预处理的STAR索引路径
        # star_index_path 已在前面从 genome_paths 中获取

        # 4) 构造 sample_inputs（直接使用 FastP 返回的结构）
        sample_inputs: List[Dict[str, Any]] = []
        for index, fastp_sample_output in enumerate(fastp_per_sample_outputs):
            sample_id = fastp_sample_output.get("sample_id", f"sample_{index + 1}")
            trimmed_read1_path = fastp_sample_output.get("trimmed_single") or fastp_sample_output.get("trimmed_r1")
            trimmed_read2_path = fastp_sample_output.get("trimmed_r2")
            if not trimmed_read1_path:
                continue
            is_paired_sample = fastp_sample_output.get("paired_end")
            if is_paired_sample is None:
                is_paired_sample = bool(trimmed_read2_path)
            else:
                is_paired_sample = bool(is_paired_sample)
            sample_input_entry = {
                "sample_id": sample_id,
                "is_paired": is_paired_sample,
                "read1": trimmed_read1_path,
            }
            if is_paired_sample and trimmed_read2_path:
                sample_input_entry["read2"] = trimmed_read2_path
            elif not is_paired_sample and trimmed_read2_path:
                # 兜底：fastp 结果声明为单端但仍给出了 read2，保留文件以便排查
                sample_input_entry["read2"] = trimmed_read2_path
            sample_inputs.append(sample_input_entry)
        if not sample_inputs:
            return {"success": False, "error": "未从FastP结果构造到任何样本输入"}

        # 5) 组装 Nextflow 参数
        cleaned_params: Dict[str, Any] = {}
        for parameter_name, parameter_value in (star_params or {}).items():
            if parameter_value is None or parameter_name in {"star_cpus", "outFileNamePrefix"}:
                continue
            cleaned_params[parameter_name.lstrip('-')] = parameter_value

        nextflow_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            "star_index": str(star_index_path),
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **cleaned_params,
        }

        # 资源配置：直接内联通过 -params 传入
        from .utils_tools import build_stage_resources_map
        resource_config_map = resource_config or {}
        resources_map = build_stage_resources_map(resource_config_map, ["star"])

        # 参数版本化
        try:
            from .utils_tools import write_params_file
            versioned_params_file = write_params_file(
                "star",
                nextflow_params,
                results_dir=str(results_dir)
            )

            if not versioned_params_file or not versioned_params_file.exists():
                raise Exception("版本化参数文件创建失败")

            logger.info(f"STAR参数版本化文件: {versioned_params_file}")
            params_file = versioned_params_file

        except Exception as exception:
            return {"success": False, "error": f"STAR参数版本化失败: {exception}"}

        # 6) 定位并执行 Nextflow
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "star.nf"
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 star.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "star.nf")]
            }

        logger.info(f"执行STAR比对 - 参数文件: {params_file}")
        logger.info(f"STAR索引: {nextflow_params['star_index']}")
        command = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        # 通过 -params 内联注入资源配置
        try:
            inline_params = json.dumps({"resources": resources_map}, ensure_ascii=False)
            command.extend(["-params", inline_params])
        except Exception as e:
            logger.warning(f"构建STAR内联资源参数失败，将不注入资源: {e}")
        execution_result = subprocess.run(command, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 star.nf publishDir 对齐）
        star_output_dir = results_dir / "star"
        per_sample_outputs: List[Dict[str, Any]] = []
        for sample_input in sample_inputs:
            sample_id = sample_input["sample_id"]
            sample_dir = star_output_dir / sample_id
            sample_result_entry = {
                "sample_id": sample_id,
                "aligned_bam": str(sample_dir / f"{sample_id}.Aligned.sortedByCoord.out.bam"),
                "log_final": str(sample_dir / f"{sample_id}.Log.final.out"),
                "log_out": str(sample_dir / f"{sample_id}.Log.out"),
                "log_progress": str(sample_dir / f"{sample_id}.Log.progress.out"),
                "splice_junctions": str(sample_dir / f"{sample_id}.SJ.out.tab"),
            }
            quant_mode = str(nextflow_params.get("quantMode", ""))
            if "TranscriptomeSAM" in quant_mode:
                sample_result_entry["transcriptome_bam"] = str(sample_dir / f"{sample_id}.Aligned.toTranscriptome.out.bam")
            if "GeneCounts" in quant_mode:
                sample_result_entry["gene_counts"] = str(sample_dir / f"{sample_id}.ReadsPerGene.out.tab")
            per_sample_outputs.append(sample_result_entry)

        results = {
            "success": execution_result.returncode == 0,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": len(sample_inputs),
            "per_sample_outputs": per_sample_outputs,
        }
        if get_tools_config().settings.debug_mode:
            results.update({"stderr": execution_result.stderr, "cmd": " ".join(command)})
        if results["success"]:
            logger.info(f"STAR完成:samples={len(sample_inputs)} results={results_dir}")
        else:
            logger.warning(f"STAR失败:rc={execution_result.returncode} stderr={(execution_result.stderr or '')[:400]}")
            # 失败时也添加调试信息
            results.update({
                "return_code": execution_result.returncode,
                "stderr": execution_result.stderr[:1000] if execution_result.stderr else "",
                "cmd": " ".join(command),
                "work_dir": str(work_dir)
            })
        return results

    except Exception as exception:
        logger.error(f"STAR异常：{exception}")
        return {"success": False, "error": f"执行STAR比对失败: {str(exception)}"}


def _parse_star_log_final(log_file_path: Path) -> Dict[str, Any]:
    """解析单个STAR Log.final.out文件，提取所有参数

    Args:
        log_file_path: Log.final.out文件路径

    Returns:
        Dict: 解析后的所有STAR参数
    """
    try:
        with open(log_file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        metrics = {}

        # 按行解析，用"|"分割每行
        for line in content.split('\n'):
            line = line.strip()
            if '|' in line:
                # 按"|"分割为两部分
                parts = line.split('|', 1)  # 只分割第一个"|"
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()

                    # 跳过空值
                    if not key or not value:
                        continue

                    # 解析数值
                    parsed_value = _parse_value(value)
                    metrics[key] = parsed_value

        return metrics

    except Exception as e:
        logger.error(f"解析STAR日志文件失败 {log_file_path}: {e}")
        return {}


def _parse_value(value_str: str) -> Any:
    """解析STAR日志中的值，自动识别类型

    Args:
        value_str: 原始字符串值

    Returns:
        解析后的值（int、float或str）
    """
    value_str = value_str.strip()

    # 处理百分比
    if value_str.endswith('%'):
        try:
            return float(value_str[:-1])
        except:
            return value_str

    # 处理数字（可能包含逗号）
    clean_value = value_str.replace(',', '')

    # 尝试解析为整数
    try:
        return int(clean_value)
    except:
        pass

    # 尝试解析为浮点数
    try:
        return float(clean_value)
    except:
        pass

    # 返回原始字符串
    return value_str


def _extract_key_metrics(all_metrics: Dict[str, Any]) -> Dict[str, Any]:
    """从所有解析的参数中提取关键指标

    Args:
        all_metrics: 解析后的所有STAR参数

    Returns:
        Dict: 标准化的关键指标
    """
    # 定义关键指标的映射关系
    key_mappings = {
        # 输入reads
        "Number of input reads": "input_reads",

        # 唯一比对
        "Uniquely mapped reads number": "uniquely_mapped",
        "Uniquely mapped reads %": "uniquely_mapped_percentage",

        # 错配率
        "Mismatch rate per base, %": "mismatch_rate",

        # 多重比对
        "Number of reads mapped to multiple loci": "multimapped",
        "% of reads mapped to multiple loci": "multimapped_percentage",

        # 未比对reads
        "Number of reads unmapped: too many mismatches": "unmapped_mismatches",
        "% of reads unmapped: too many mismatches": "unmapped_mismatches_percentage",
        "Number of reads unmapped: too short": "unmapped_too_short",
        "% of reads unmapped: too short": "unmapped_too_short_percentage",
        "Number of reads unmapped: other": "unmapped_other",
        "% of reads unmapped: other": "unmapped_other_percentage",

        # 剪接相关
        "Number of splices: Total": "total_splices",
        "Number of splices: Annotated (sjdb)": "annotated_splices",
        "Number of splices: GT/AG": "gt_ag_splices",
        "Number of splices: GC/AG": "gc_ag_splices",
        "Number of splices: AT/AC": "at_ac_splices",
        "Number of splices: Non-canonical": "non_canonical_splices",

        # 多重map reads
        "Number of reads mapped to too many loci": "reads_too_many_loci",
        "% of reads mapped to too many loci": "reads_too_many_loci_percentage",

        # 删除/插入
        "Deletion rate per base": "deletion_rate",
        "Deletion average length": "deletion_avg_length",
        "Insertion rate per base": "insertion_rate",
        "Insertion average length": "insertion_avg_length"
    }

    # 提取关键指标
    key_metrics = {}
    for original_key, standard_key in key_mappings.items():
        if original_key in all_metrics:
            key_metrics[standard_key] = all_metrics[original_key]

    # 计算总的未比对数量和百分比
    if "input_reads" in key_metrics:
        total_unmapped = 0
        unmapped_keys = ["unmapped_mismatches", "unmapped_too_short", "unmapped_other"]

        for key in unmapped_keys:
            if key in key_metrics:
                total_unmapped += key_metrics[key]

        key_metrics["unmapped"] = total_unmapped

        if key_metrics["input_reads"] > 0:
            key_metrics["unmapped_percentage"] = round(
                (total_unmapped / key_metrics["input_reads"]) * 100, 2
            )

    return key_metrics


@tool
def parse_star_metrics(results_directory: str) -> Dict[str, Any]:
    """解析STAR比对结果文件，深入提取所有统计指标

    Args:
        results_directory: STAR结果目录路径

    Returns:
        Dict: 包含完整统计指标的数据
        {
            "success": bool,
            "total_samples": int,
            "overall_statistics": Dict,  # 汇总统计
            "sample_metrics": List[Dict],  # 每个样本的详细指标
            "raw_metrics": Dict,  # 原始解析的所有参数
            "results_directory": str
        }
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {
                "success": False,
                "error": f"STAR结果目录不存在: {results_directory}"
            }

        # 查找STAR输出目录
        star_dir = results_path / "star"
        if not star_dir.exists():
            return {
                "success": False,
                "error": f"STAR输出目录不存在: {star_dir}"
            }

        sample_metrics = []
        all_raw_metrics = {}
        overall_stats = {
            "total_input_reads": 0,
            "total_uniquely_mapped": 0,
            "total_multimapped": 0,
            "total_unmapped": 0,
            "overall_uniquely_mapped_percentage": 0.0,
            "overall_multimapped_percentage": 0.0,
            "overall_mapped_percentage": 0.0
        }

        # 扫描每个样本的STAR输出
        sample_dirs = [d for d in star_dir.iterdir() if d.is_dir()]

        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            log_final_path = sample_dir / "Log.final.out"

            if not log_final_path.exists():
                logger.warning(f"未找到STAR日志文件: {log_final_path}")
                continue

            # 解析完整的日志文件
            raw_metrics = _parse_star_log_final(log_final_path)
            if not raw_metrics:
                continue

            # 提取关键指标
            key_metrics = _extract_key_metrics(raw_metrics)
            key_metrics["sample_id"] = sample_id

            sample_metrics.append(key_metrics)
            all_raw_metrics[sample_id] = raw_metrics

            # 累计到总体统计
            if "input_reads" in key_metrics:
                overall_stats["total_input_reads"] += key_metrics["input_reads"]
            if "uniquely_mapped" in key_metrics:
                overall_stats["total_uniquely_mapped"] += key_metrics["uniquely_mapped"]
            if "multimapped" in key_metrics:
                overall_stats["total_multimapped"] += key_metrics["multimapped"]
            if "unmapped" in key_metrics:
                overall_stats["total_unmapped"] += key_metrics["unmapped"]

        # 计算总体百分比
        if overall_stats["total_input_reads"] > 0:
            total_reads = overall_stats["total_input_reads"]
            overall_stats["overall_uniquely_mapped_percentage"] = round(
                (overall_stats["total_uniquely_mapped"] / total_reads) * 100, 2
            )
            overall_stats["overall_multimapped_percentage"] = round(
                (overall_stats["total_multimapped"] / total_reads) * 100, 2
            )
            overall_stats["overall_mapped_percentage"] = round(
                ((overall_stats["total_uniquely_mapped"] + overall_stats["total_multimapped"]) / total_reads) * 100, 2
            )

        return {
            "success": True,
            "total_samples": len(sample_metrics),
            "overall_statistics": overall_stats,
            "sample_metrics": sample_metrics,
            "raw_metrics": all_raw_metrics,  # 新增：包含所有原始解析参数
            "results_directory": results_directory
        }

    except Exception as e:
        logger.error(f"解析STAR指标失败: {e}")
        return {
            "success": False,
            "error": f"解析STAR指标失败: {str(e)}"
        }
