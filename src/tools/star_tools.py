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
    genome_id: str
) -> Dict[str, Any]:
    """执行 STAR 比对（精简版）

    Args:
        star_params: STAR执行参数
        fastp_results: FastP质控结果
        genome_id: 基因组ID（如"dm6", "hg38"等）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - 基于genome_id从genomes.json获取配置并构建STAR索引路径
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 star.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 读取genomes.json获取基因组配置
        genomes_config_path = tools_config.genomes_config_path
        if not genomes_config_path.exists():
            return {"success": False, "error": f"基因组配置文件不存在: {genomes_config_path}"}

        try:
            with open(genomes_config_path, 'r', encoding='utf-8') as f:
                genomes_config = json.load(f)
        except Exception as exception:
            return {"success": False, "error": f"读取基因组配置失败: {exception}"}

        if genome_id not in genomes_config:
            return {"success": False, "error": f"基因组配置中不存在: {genome_id}"}

        genome_info = genomes_config[genome_id]

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

        # 3) 直接构建STAR索引路径（detect节点已验证索引存在性）
        species = genome_info.get("species")
        version = genome_info.get("version")
        if not species or not version:
            return {"success": False, "error": f"基因组配置缺少species或version字段: {genome_id}"}

        star_index_path = tools_config.settings.data_dir / "genomes" / species / version / "star_index"

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


@tool
def parse_star_metrics(results_directory: str) -> Dict[str, Any]:
    """解析STAR比对结果文件，返回原始日志文件内容
    
    Args:
        results_directory: STAR结果目录路径
    
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
                "error": f"STAR结果目录不存在: {results_directory}"
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
        
        sample_logs = []
        
        # 遍历样本目录，收集Log.final.out文件内容
        for sample_dir in star_dir.iterdir():
            if not sample_dir.is_dir():
                continue
                
            sample_id = sample_dir.name
            log_final_file = sample_dir / f"{sample_id}.Log.final.out"
            
            if not log_final_file.exists():
                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(log_final_file),
                    "error": "Log.final.out文件不存在"
                })
                continue
            
            try:
                # 读取Log.final.out原始内容
                with open(log_final_file, 'r', encoding='utf-8') as f:
                    log_content = f.read()
                
                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(log_final_file),
                    "content": log_content
                })
                
            except Exception as exception:
                sample_logs.append({
                    "sample_id": sample_id,
                    "log_file": str(log_final_file),
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
            "sample_logs": sample_logs
        }
        
        try:
            logger.info(f"STAR日志收集完成: {successful_samples}个样本，目录={results_directory}")
        except Exception:
            pass
        
        return result
        
    except Exception as exception:
        try:
            logger.error(f"收集STAR日志失败：{exception}")
        except Exception:
            pass
        return {
            "success": False,
            "error": f"收集STAR日志失败: {str(exception)}"
        }
