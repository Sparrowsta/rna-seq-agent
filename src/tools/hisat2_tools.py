"""
RNA-seq智能分析助手 - HISAT2工具集

包含：
- run_nextflow_hisat2: 执行HISAT2比对
- parse_hisat2_metrics: 解析HISAT2结果指标
"""

import json
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.hisat2")


@tool
def run_nextflow_hisat2(
    hisat2_params: Dict[str, Any],
    fastp_results: Dict[str, Any],
    genome_id: str,
    resource_config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """执行 HISAT2 比对

    Args:
        hisat2_params: HISAT2执行参数
        fastp_results: FastP质控结果
        genome_id: 基因组ID，用于获取对应的基因组路径信息

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - 通过 genome_id 从配置文件中获取HISAT2索引路径
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs
    - per_sample_outputs 路径与 hisat2.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 通过 genome_id 获取基因组路径信息
        from .utils_tools import get_genome_paths_from_id
        genome_paths = get_genome_paths_from_id(genome_id)

        if not genome_paths:
            return {
                "success": False,
                "error": f"无法获取基因组路径信息: {genome_id}"
            }

        hisat2_index_path = genome_paths.get("hisat2_index_path", "")
        if not hisat2_index_path:
            return {
                "success": False,
                "error": f"基因组配置缺少HISAT2索引路径: {genome_id}"
            }

        # 2) 根据hisat2_index_path推断hisat2_index_prefix
        try:
            import re
            idx_path = Path(hisat2_index_path)
            if not idx_path.is_absolute():
                idx_path = tools_config.settings.data_dir / idx_path
            if idx_path.is_dir():
                candidates = sorted(idx_path.glob('*.1.ht2'))
                if not candidates:
                    candidates = sorted(idx_path.glob('*.ht2'))
                if candidates:
                    first_name = candidates[0].name
                    prefix_base = re.sub(r'\.(?:[1-8])\.ht2$', '', first_name)
                    hisat2_index_prefix = str(idx_path / prefix_base)
                else:
                    hisat2_index_prefix = str(idx_path / 'genome')
            else:
                hisat2_index_prefix = str(idx_path)
        except Exception as _e:
            logger.warning(f"HISAT2 索引前缀推断失败，使用原始配置: {hisat2_index_path} (error={_e})")
            hisat2_index_prefix = str(hisat2_index_path)

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行HISAT2比对"}

        fastp_results_dir = fastp_results.get("results_dir")
        if not fastp_results_dir:
            return {"success": False, "error": "FastP结果缺少results_dir"}

        fastp_per_sample_outputs = fastp_results.get("per_sample_outputs") or []
        if not fastp_per_sample_outputs:
            return {"success": False, "error": "FastP结果缺少per_sample_outputs"}

        # 2) 运行根目录与工作目录
        results_dir = Path(fastp_results_dir)
        run_id = results_dir.name
        work_dir = tools_config.settings.data_dir / "work" / f"hisat2_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"HISAT2启动:samples={len(fastp_per_sample_outputs)} results={results_dir} work={work_dir}")

        # 3) 直接使用预处理的HISAT2索引路径
        # hisat2_index_path 已在前面从 genome_paths 中获取
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

        # 5) 资源与 Nextflow 参数
        # 先归一化并构建资源映射，再写入参数文件，避免引用顺序问题
        from .utils_tools import normalize_resources
        normalized_hisat2 = normalize_resources("hisat2", {"hisat2": resource_config or {}})
        resources_map: Dict[str, Dict[str, Any]] = {"hisat2": normalized_hisat2} if normalized_hisat2 else {}

        cleaned_params: Dict[str, Any] = {}
        for parameter_name, parameter_value in (hisat2_params or {}).items():
            if parameter_value is None or parameter_name in {"hisat2_cpus", "threads", "p"}:
                continue
            cleaned_params[parameter_name.lstrip('-')] = parameter_value

        nextflow_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            # 传递给 Nextflow 的应为索引前缀的绝对路径
            "hisat2_index": str(hisat2_index_prefix),
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            # 将资源配置并入 params-file，避免对 -params 的版本依赖
            "resources": resources_map,
            **cleaned_params,
        }

        # 参数版本化
        try:
            from .utils_tools import write_params_file
            versioned_params_file = write_params_file(
                "hisat2",
                nextflow_params,
                results_dir=str(results_dir)
            )

            if not versioned_params_file or not versioned_params_file.exists():
                raise Exception("版本化参数文件创建失败")

            logger.info(f"HISAT2参数版本化文件: {versioned_params_file}")
            params_file = versioned_params_file

        except Exception as exception:
            return {"success": False, "error": f"HISAT2参数版本化失败: {exception}"}

        # 6) 定位并执行 Nextflow
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "hisat2.nf"
        if not nextflow_script.exists():
            return {
                "success": False, 
                "error": "未找到 hisat2.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "hisat2.nf")]
            }

        logger.info(f"执行HISAT2比对 - 参数文件: {params_file}")
        logger.info(f"HISAT2索引前缀: {nextflow_params.get('hisat2_index')}")
        command = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        # 统一通过 params-file 传递资源，避免内联 -params
        execution_result = subprocess.run(command, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 hisat2.nf publishDir 对齐）
        hisat2_output_dir = results_dir / "hisat2"
        per_sample_outputs: List[Dict[str, Any]] = []
        for sample_input in sample_inputs:
            sample_id = sample_input["sample_id"]
            sample_dir = hisat2_output_dir / sample_id
            sample_result_entry = {
                "sample_id": sample_id,
                "aligned_bam": str(sample_dir / f"{sample_id}.hisat2.bam"),
                "align_summary": str(sample_dir / f"{sample_id}.align_summary.txt"),
                "bam_index": str(sample_dir / f"{sample_id}.hisat2.bam.bai"),
            }
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
            logger.info(f"HISAT2完成:samples={len(sample_inputs)} results={results_dir}")
        else:
            logger.warning(f"HISAT2失败:rc={execution_result.returncode} stderr={(execution_result.stderr or '')[:400]}")
            # 失败时也添加调试信息
            results.update({
                "return_code": execution_result.returncode,
                "stderr": execution_result.stderr[:1000] if execution_result.stderr else "",
                "cmd": " ".join(command),
                "work_dir": str(work_dir)
            })
        return results

    except Exception as exception:
        logger.error(f"HISAT2异常：{exception}")
        return {"success": False, "error": f"执行HISAT2比对失败: {str(exception)}"}


@tool
def parse_hisat2_metrics(results_directory: str) -> Dict[str, Any]:
    """解析HISAT2比对结果文件，返回原始日志文件内容

    Args:
        results_directory: HISAT2结果目录路径

    Returns:
        Dict: 包含原始日志文件内容的数据
        {
            "success": bool,
            "total_samples": int,
            "results_directory": str,
            "sample_logs": List[Dict]  # 每个样本的原始日志内容
        }
    """
    try:
        results_path = Path(results_directory)

        if not results_path.exists():
            return {
                "success": False,
                "error": f"HISAT2结果目录不存在: {results_directory}"
            }

        # 查找HISAT2输出目录
        hisat2_dir = results_path / "hisat2"
        if not hisat2_dir.exists():
            return {
                "success": False,
                "error": f"HISAT2输出目录不存在: {hisat2_dir}"
            }

        sample_logs = []
        sample_metrics: List[Dict[str, Any]] = []

        # 遍历样本目录，收集align_summary.txt文件内容
        for sample_dir in hisat2_dir.iterdir():
            if not sample_dir.is_dir():
                continue

            sample_id = sample_dir.name
            summary_file = sample_dir / f"{sample_id}.align_summary.txt"

            if not summary_file.exists():
                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(summary_file),
                    "error": "align_summary.txt文件不存在"
                })
                continue

            try:
                # 读取align_summary.txt原始内容
                with open(summary_file, 'r', encoding='utf-8') as f:
                    log_content = f.read()

                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(summary_file),
                    "content": log_content
                })

                # 解析关键指标：overall alignment rate、总reads数等
                import re
                overall_rate = None
                total_reads = None

                # 匹配例如："123456 reads; of these:" 或 "123,456 reads; of these:"
                m_total = re.search(r"([\d,]+)\s+reads;\s+of\s+these:", log_content, re.IGNORECASE)
                if m_total:
                    try:
                        total_reads = int(m_total.group(1).replace(',', ''))
                    except Exception:
                        total_reads = None

                # 匹配例如："95.67% overall alignment rate"
                m_rate = re.search(r"([\d.]+)%\s+overall\s+alignment\s+rate", log_content, re.IGNORECASE)
                if m_rate:
                    try:
                        overall_rate = float(m_rate.group(1))
                    except Exception:
                        overall_rate = None

                metric_entry = {"sample_id": sample_id}
                if overall_rate is not None:
                    metric_entry["overall_alignment_rate"] = overall_rate
                if total_reads is not None:
                    metric_entry["total_reads"] = total_reads

                # 仅在至少有一项可用指标时加入
                if len(metric_entry) > 1:
                    sample_metrics.append(metric_entry)

            except Exception as exception:
                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(summary_file),
                    "error": f"读取失败: {str(exception)}"
                })

        if not sample_logs:
            return {
                "success": False,
                "error": "未找到任何样本日志文件"
            }

        # 计算成功读取的样本数
        successful_samples = len([log for log in sample_logs if "error" not in log])

        result = {
            "success": True,
            "total_samples": successful_samples,
            "results_directory": results_directory,
            "sample_logs": sample_logs,
            # 新增：为对齐函数提供可消费的样本级指标
            "sample_metrics": sample_metrics,
        }

        try:
            logger.info(f"HISAT2日志收集完成: {successful_samples}个样本，目录={results_directory}")
        except Exception:
            pass

        return result

    except Exception as exception:
        try:
            logger.error(f"收集HISAT2日志失败：{exception}")
        except Exception:
            pass
        return {
            "success": False,
            "error": f"收集HISAT2日志失败: {str(exception)}"
        }
