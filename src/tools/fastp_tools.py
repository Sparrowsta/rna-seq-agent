"""
RNA-seq智能分析助手 - FastP质控工具

包含：
- run_nextflow_fastp: 执行FastP质控流程
- parse_fastp_results: 解析FastP结果
"""

import json
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.fastp")


@tool
def run_nextflow_fastp(fastp_params: Dict[str, Any], sample_info: Dict[str, Any]) -> Dict[str, Any]:
    """执行Nextflow FastP质控流程
    
    Args:
        fastp_params: FastP参数字典，例如 {"qualified_quality_phred": 25, "length_required": 50}
        sample_info: 样本信息，包含sample_groups、state_info等
    
    Returns:
        执行结果字典，包含状态、输出路径、执行日志等
    """
    try:
        config = get_tools_config()
        
        # 校验必需参数
        if not fastp_params:
            return {
                "success": False,
                "error": "FastP参数不能为空",
                "execution_time": 0
            }
        
        if not sample_info.get("sample_groups"):
            return {
                "success": False, 
                "error": "样本信息缺失",
                "execution_time": 0
            }
        
        # 记录开始时间
        start_time = time.time()
        
        # 统一数据根目录来源：始终以 Settings().data_dir 为准，不从 sample_info 读取
        base_data_path = str(config.settings.data_dir)

        # 规范化 paired_end：优先使用显式传入，其次从 sample_groups 是否存在 read2 推断；默认 False
        paired_end_value = sample_info.get("paired_end")
        if paired_end_value is None:
            groups = sample_info.get("sample_groups", [])
            paired_end_value = any((isinstance(s, dict) and s.get("read2")) for s in groups)
        
        # 结果目录（运行根目录）：优先使用 sample_info 提供；否则按时间戳生成
        if sample_info and sample_info.get("results_dir"):
            results_dir = Path(sample_info["results_dir"])
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            results_dir = config.settings.data_dir / "results" / f"fastp_{timestamp}"

        # 工作目录：统一到 /data/work
        run_id = results_dir.name
        work_dir = config.settings.data_dir / 'work' / f"fastp_{run_id}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 确保结果目录存在，避免 publishDir 目标不存在造成的发布失败
        try:
            results_dir.mkdir(parents=True, exist_ok=True)
            (results_dir / "fastp").mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.warning(f"创建结果目录失败: {e}")
        
        # 准备Nextflow配置（用于 -params-file）
        nextflow_config = {
            "work_dir": str(work_dir),
            "results_dir": str(results_dir),
            "data_dir": base_data_path,
            "sample_groups": sample_info["sample_groups"],
            "paired_end": bool(paired_end_value),
            **fastp_params
        }

        # 写入参数版本化文件
        versioned_params_file = None
        try:
            from .utils_tools import write_params_file
            # 直接传递results_dir，不需要创建临时state
            versioned_params_file = write_params_file(
                "fastp",
                nextflow_config,
                results_dir=str(results_dir)
            )
            logger.info(f"FastP实际运行参数版本化文件已写入: {versioned_params_file}")

        except Exception as e:
            logger.warning(f"FastP参数版本化写入失败 (不影响执行): {e}")

        # 使用版本化文件作为参数文件，如果失败则降级到临时文件
        if versioned_params_file and versioned_params_file.exists():
            params_file = versioned_params_file
            logger.info(f"使用版本化参数文件运行FastP: {params_file}")
        else:
            # 降级方案：生成临时参数文件
            logger.warning("版本化参数文件不可用，使用临时参数文件")
            params_file = work_dir / "nextflow.config.json"
            with open(params_file, 'w', encoding='utf-8') as f:
                json.dump(nextflow_config, f, indent=2, ensure_ascii=False)
        
        # 寻找Nextflow脚本
        nextflow_script = config.settings.nextflow_scripts_dir / "fastp.nf"
        
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 fastp.nf",
                "searched": [str(config.settings.nextflow_scripts_dir / "fastp.nf")],
                "execution_time": 0
            }
        
        # 执行Nextflow
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
            "-resume"
        ]
        
        logger.info(f"执行FastP: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800,  # 30分钟超时
                cwd=str(work_dir)
            )
            
            execution_time = time.time() - start_time
            
            if result.returncode == 0:
                # 解析输出文件
                per_sample_outputs = []
                for sample in sample_info["sample_groups"]:
                    sample_id = sample["sample_id"]
                    sample_dir = results_dir / "fastp" / sample_id
                    
                    output_info = {
                        "sample_id": sample_id,
                        "html": str(sample_dir / f"{sample_id}.fastp.html"),
                        "json": str(sample_dir / f"{sample_id}.fastp.json")
                    }
                    
                    if paired_end_value:
                        output_info.update({
                            "trimmed_r1": str(sample_dir / f"{sample_id}_1.trimmed.fastq.gz"),
                            "trimmed_r2": str(sample_dir / f"{sample_id}_2.trimmed.fastq.gz")
                        })
                    else:
                        output_info["trimmed_single"] = str(sample_dir / f"{sample_id}.single.trimmed.fastq.gz")
                    
                    per_sample_outputs.append(output_info)
                
                return {
                    "success": True,
                    "results_dir": str(results_dir),
                    "per_sample_outputs": per_sample_outputs,
                    "execution_time": execution_time,
                    "params_file": str(params_file),
                    "message": "FastP执行成功"
                }
            else:
                return {
                    "success": False,
                    "error": f"FastP执行失败，返回码: {result.returncode}",
                    "stderr": result.stderr,
                    "execution_time": execution_time
                }
        
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": "FastP执行超时（30分钟）",
                "execution_time": time.time() - start_time
            }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"FastP执行异常: {str(e)}",
            "execution_time": 0
        }


@tool
def parse_fastp_results(results_directory: str) -> Dict[str, Any]:
    """解析FastP结果文件，提取基本质量统计和特殊处理的复杂数据供LLM分析

    Args:
        results_directory: FastP结果目录路径

    Returns:
        解析的质量指标字典，包含：
        - 基本前后质量统计（reads数量、质量分数、GC含量等）
        - Top10 k-mer和过度表达序列
        - 原始复杂数据直接传入LLM分析
    """
    try:
        results_dir = Path(results_directory)
        if not results_dir.exists():
            return {
                "success": False,
                "error": f"结果目录不存在: {results_directory}"
            }

        # 查找所有FastP JSON报告文件
        json_files = list(results_dir.rglob("*.fastp.json"))

        if not json_files:
            logger.warning(f"未找到FastP JSON报告文件：{results_directory}")
            return {
                "success": False,
                "error": "未找到FastP JSON报告文件"
            }

        sample_metrics = []
        overall_stats = {
            "total_reads_before": 0,
            "total_reads_after": 0,
            "total_bases_before": 0,
            "total_bases_after": 0,
            "q20_rates": [],
            "q30_rates": [],
            "gc_contents": []
        }

        # 解析每个样本的JSON报告
        for json_file in json_files:
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    fastp_data = json.load(f)

                # 提取样本ID
                sample_id = json_file.stem.replace('.fastp', '')

                # 基本统计信息
                summary = fastp_data.get("summary", {})
                before_filtering = summary.get("before_filtering", {})
                after_filtering = summary.get("after_filtering", {})

                # 基本质量统计（前后对比）
                reads_before = before_filtering.get("total_reads", 0)
                reads_after = after_filtering.get("total_reads", 0)
                bases_before = before_filtering.get("total_bases", 0)
                bases_after = after_filtering.get("total_bases", 0)

                # 质量分数统计
                q20_before = before_filtering.get("q20_bases", 0)
                q30_before = before_filtering.get("q30_bases", 0)
                q20_after = after_filtering.get("q20_bases", 0)
                q30_after = after_filtering.get("q30_bases", 0)

                # 计算质量率
                q20_rate_before = (q20_before / bases_before * 100) if bases_before > 0 else 0
                q30_rate_before = (q30_before / bases_before * 100) if bases_before > 0 else 0
                q20_rate_after = (q20_after / bases_after * 100) if bases_after > 0 else 0
                q30_rate_after = (q30_after / bases_after * 100) if bases_after > 0 else 0

                # GC含量
                gc_before = before_filtering.get("gc_content", 0)
                gc_after = after_filtering.get("gc_content", 0)

                # 读长统计
                read1_mean_length_before = before_filtering.get("read1_mean_length", 0)
                read1_mean_length_after = after_filtering.get("read1_mean_length", 0)

                # 特殊处理：计算前后变化最大的Top10 k-mer
                top_changed_kmers = []
                kmer_data_before = fastp_data.get("read1_before_filtering", {}).get("kmer_count", {})
                kmer_data_after = fastp_data.get("read1_after_filtering", {}).get("kmer_count", {})

                if kmer_data_before and kmer_data_after:
                    kmer_changes = []
                    # 计算所有k-mer的前后变化
                    for seq in kmer_data_before.keys():
                        count_before = kmer_data_before[seq]
                        count_after = kmer_data_after.get(seq, 0)  # 如果after中没有，说明被完全去除
                        absolute_change = abs(count_before - count_after)
                        change_rate = ((count_before - count_after) / count_before * 100) if count_before > 0 else 0

                        kmer_changes.append({
                            "sequence": seq,
                            "count_before": count_before,
                            "count_after": count_after,
                            "absolute_change": absolute_change,
                            "change_rate": round(change_rate, 2)
                        })

                    # 按绝对变化量降序排序，取前10个
                    kmer_changes.sort(key=lambda x: x["absolute_change"], reverse=True)
                    top_changed_kmers = kmer_changes[:10]

                # 特殊处理：计算前后变化最大的Top10 overrepresented_sequences
                top_changed_overrepresented = []
                overrep_data_before = fastp_data.get("read1_before_filtering", {}).get("overrepresented_sequences", {})
                overrep_data_after = fastp_data.get("read1_after_filtering", {}).get("overrepresented_sequences", {})

                if overrep_data_before and overrep_data_after:
                    overrep_changes = []
                    # 计算所有过度表达序列的前后变化
                    for seq in overrep_data_before.keys():
                        count_before = overrep_data_before[seq]
                        count_after = overrep_data_after.get(seq, 0)  # 如果after中没有，说明被完全去除
                        absolute_change = abs(count_before - count_after)
                        change_rate = ((count_before - count_after) / count_before * 100) if count_before > 0 else 0

                        overrep_changes.append({
                            "sequence": seq,
                            "count_before": count_before,
                            "count_after": count_after,
                            "absolute_change": absolute_change,
                            "change_rate": round(change_rate, 2)
                        })

                    # 按绝对变化量降序排序，取前10个
                    overrep_changes.sort(key=lambda x: x["absolute_change"], reverse=True)
                    top_changed_overrepresented = overrep_changes[:10]

                sample_metric = {
                    "sample_id": sample_id,
                    # 基本前后统计
                    "reads_before": reads_before,
                    "reads_after": reads_after,
                    "reads_passed_rate": round((reads_after / reads_before * 100) if reads_before > 0 else 0, 2),
                    "bases_before": bases_before,
                    "bases_after": bases_after,
                    "bases_passed_rate": round((bases_after / bases_before * 100) if bases_before > 0 else 0, 2),
                    "q20_rate_before": round(q20_rate_before, 2),
                    "q20_rate_after": round(q20_rate_after, 2),
                    "q30_rate_before": round(q30_rate_before, 2),
                    "q30_rate_after": round(q30_rate_after, 2),
                    "gc_content_before": round(gc_before, 2),
                    "gc_content_after": round(gc_after, 2),
                    "mean_length_before": read1_mean_length_before,
                    "mean_length_after": read1_mean_length_after,

                    # 特殊处理的Top10变化数据
                    "top_changed_kmers": top_changed_kmers,
                    "top_changed_overrepresented": top_changed_overrepresented,

                    # 直接传入LLM的原始复杂数据
                    "raw_filtering_result": fastp_data.get("filtering_result", {}),
                    "raw_adapter_cutting": fastp_data.get("adapter_cutting", {}),
                    "raw_duplication": fastp_data.get("duplication", {}),

                    "json_file": str(json_file)
                }

                sample_metrics.append(sample_metric)

                # 累计统计
                overall_stats["total_reads_before"] += reads_before
                overall_stats["total_reads_after"] += reads_after
                overall_stats["total_bases_before"] += bases_before
                overall_stats["total_bases_after"] += bases_after
                overall_stats["q20_rates"].append(q20_rate_after)
                overall_stats["q30_rates"].append(q30_rate_after)
                overall_stats["gc_contents"].append(gc_after)

            except Exception as e:
                logger.error(f"解析FastP JSON文件失败: {json_file}, 错误: {e}")
                continue

        if not sample_metrics:
            return {
                "success": False,
                "error": "没有成功解析任何FastP结果文件"
            }

        # 计算整体统计
        total_samples = len(sample_metrics)
        overall_read_retention = (overall_stats["total_reads_after"] / overall_stats["total_reads_before"] * 100) if overall_stats["total_reads_before"] > 0 else 0
        overall_base_retention = (overall_stats["total_bases_after"] / overall_stats["total_bases_before"] * 100) if overall_stats["total_bases_before"] > 0 else 0

        avg_q20_rate = sum(overall_stats["q20_rates"]) / len(overall_stats["q20_rates"]) if overall_stats["q20_rates"] else 0
        avg_q30_rate = sum(overall_stats["q30_rates"]) / len(overall_stats["q30_rates"]) if overall_stats["q30_rates"] else 0
        avg_gc_content = sum(overall_stats["gc_contents"]) / len(overall_stats["gc_contents"]) if overall_stats["gc_contents"] else 0

        result = {
            "success": True,
            "sample_count": total_samples,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_reads_before": overall_stats["total_reads_before"],
                "total_reads_after": overall_stats["total_reads_after"],
                "read_retention_rate": round(overall_read_retention, 2),
                "total_bases_before": overall_stats["total_bases_before"],
                "total_bases_after": overall_stats["total_bases_after"],
                "base_retention_rate": round(overall_base_retention, 2),
                "average_q20_rate": round(avg_q20_rate, 2),
                "average_q30_rate": round(avg_q30_rate, 2),
                "average_gc_content": round(avg_gc_content, 2)
            },
            "results_directory": results_directory
        }

        logger.info(f"FastP结果解析完成: {total_samples}个样本, 平均Q30={avg_q30_rate:.1f}%, 数据保留率={overall_read_retention:.1f}%")
        return result

    except Exception as e:
        logger.error(f"解析FastP结果异常: {e}")
        return {
            "success": False,
            "error": f"解析FastP结果失败: {str(e)}"
        }