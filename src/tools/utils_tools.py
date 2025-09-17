"""
RNA-seq智能分析助手 - 通用工具

包含：
- get_project_overview: 获取项目概览
- list_analysis_history: 列出分析历史
- write_params_file: 写入参数文件
- enhance_tool_result_with_debug: 增强工具结果调试信息
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
from ..state import AgentState

logger = get_logger("rna.tools.utils")


@tool
def get_project_overview() -> Dict[str, Any]:
    """获取项目整体状态概览，包括数据、基因组、系统资源和分析历史"""
    try:
        from .common_tools import scan_fastq_files, scan_system_resources, scan_genome_files
        
        # 调用各个扫描工具获取信息
        fastq_info = scan_fastq_files.invoke({})
        system_info = scan_system_resources.invoke({})
        genome_info = scan_genome_files.invoke({})
        history_info = list_analysis_history.invoke({})
        
        # 汇总概览信息
        overview = {
            "detection_status": "success",
            "timestamp": time.time(),
            "summary": {
                "fastq_files": fastq_info.get("total_files", 0),
                "fastq_samples": fastq_info.get("total_samples", 0),
                "sequencing_type": fastq_info.get("sequencing_type", "unknown"),
                "available_genomes": genome_info.get("available_genomes", 0),
                "system_memory_gb": system_info.get("memory", {}).get("total_gb", 0),
                "system_cores": system_info.get("cpu", {}).get("logical_cores", 0),
                "historical_analyses": history_info.get("total_analyses", 0)
            },
            "details": {
                "fastq_data": fastq_info,
                "system_resources": system_info,
                "genome_setup": genome_info,
                "analysis_history": history_info
            }
        }
        
        logger.info(f"项目概览: {fastq_info.get('total_samples', 0)}个样本, {genome_info.get('available_genomes', 0)}个基因组可用")
        return overview
        
    except Exception as e:
        logger.error(f"获取项目概览失败: {e}")
        return {
            "detection_status": "failed",
            "error": str(e),
            "timestamp": time.time()
        }


@tool
def list_analysis_history() -> Dict[str, Any]:
    """获取历史分析记录，返回分析时间、配置和结果信息"""
    try:
        config = get_tools_config()
        results_base = config.settings.data_dir / "results"
        
        if not results_base.exists():
            return {
                "detection_status": "no_history",
                "total_analyses": 0,
                "analyses": [],
                "message": "未找到历史分析记录"
            }
        
        analyses = []
        
        # 扫描results目录下的分析记录
        for item in results_base.iterdir():
            if not item.is_dir():
                continue
                
            # 尝试解析时间戳目录名
            try:
                timestamp = datetime.strptime(item.name, "%Y%m%d_%H%M%S")
            except ValueError:
                continue
            
            # 检查分析结果内容
            try:
                # 检查归档内容
                analysis_info = {
                    "timestamp_str": item.name,
                    "timestamp": timestamp.timestamp(),
                    "path": str(item),
                    "files": {}
                }
                
                # 检查各阶段结果
                if (item / "fastp").exists():
                    analysis_info["files"]["fastp"] = True
                if (item / "star").exists():
                    analysis_info["files"]["star"] = True
                if (item / "featurecounts").exists():
                    analysis_info["files"]["featurecounts"] = True
                if (item / "reports").exists():
                    analysis_info["files"]["reports"] = True
                
                # 计算完整度
                total_stages = 4  # fastp, star, featurecounts, reports
                completed_stages = sum(1 for v in analysis_info["files"].values() if v)
                analysis_info["completion_rate"] = (completed_stages / total_stages) * 100
                
                analyses.append(analysis_info)
                
            except Exception:
                # 忽略无法读取的目录
                continue
        
        # 按时间排序（最新在前）
        analyses.sort(key=lambda x: x["timestamp"], reverse=True)
        
        result = {
            "detection_status": "success",
            "total_analyses": len(analyses),
            "analyses": analyses[:10],  # 只返回最近10个
            "results_directory": str(results_base)
        }
        
        # 添加最新分析摘要
        if analyses:
            latest = analyses[0]
            result["latest_analysis"] = {
                "timestamp_str": latest["timestamp_str"],
                "completion_rate": latest["completion_rate"],
                "path": latest["path"]
            }
        
        # 日志
        try:
            latest = result["latest_analysis"]["timestamp_str"] if result.get("latest_analysis") else None
            logger.info(f"历史分析：total={result['total_analyses']} latest={latest}")
        except Exception:
            pass
        return result
        
    except Exception as e:
        logger.error(f"获取分析历史失败: {e}")
        return {
            "detection_status": "failed",
            "error": str(e),
            "total_analyses": 0,
            "analyses": []
        }


def write_params_file(step: str, params: dict, state: AgentState = None,
                     results_dir: str = None, metadata: dict = None) -> Path:
    """写入参数版本化文件

    Args:
        step: 执行步骤名称 (prepare/fastp/star/featurecounts等)
        params: 参数字典
        state: 当前Agent状态 (可选，向后兼容)
        results_dir: 结果目录路径 (可选，优先使用)
        metadata: 额外的元数据 (可选)

    Returns:
        写入的文件路径
    """
    try:
        config = get_tools_config()
        
        # 确定写入目录（优先使用results_dir参数，其次state.results_dir，最后默认）
        if results_dir:
            base_dir = Path(results_dir)
        elif state and hasattr(state, 'results_dir') and state.results_dir:
            base_dir = Path(state.results_dir)
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            base_dir = config.settings.data_dir / "results" / timestamp
        
        # 创建参数目录
        params_dir = base_dir / "params"
        params_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成文件名（带时间戳）
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{step}_params_{timestamp}.json"
        params_file = params_dir / filename
        
        # 准备写入的数据
        data_to_write = {
            "step": step,
            "timestamp": timestamp,
            "params": params,
            "metadata": {
                "created_by": "rna-seq-agent",
            }
        }

        # 添加state快照（如果有state）
        if state:
            data_to_write["metadata"]["agent_state_snapshot"] = {
                "status": getattr(state, 'status', ''),
                "current_step": getattr(state, 'current_step', ''),
                "execution_mode": getattr(state, 'execution_mode', '')
            }

        # 添加额外元数据（如果提供）
        if metadata:
            data_to_write["metadata"].update(metadata)
        
        # 写入文件
        with open(params_file, 'w', encoding='utf-8') as f:
            json.dump(data_to_write, f, indent=2, ensure_ascii=False)
        
        logger.info(f"参数文件已写入: {params_file}")
        return params_file
        
    except Exception as e:
        logger.error(f"写入参数文件失败: {e}")
        # 返回默认路径，避免中断流程
        fallback_path = config.settings.temp_dir / f"{step}_params_error.json"
        return fallback_path


def enhance_tool_result_with_debug(result: dict, cmd: str = "", 
                                  params_file: str = "", stdout: str = "", 
                                  stderr: str = "", execution_time: float = 0,
                                  additional_info: dict = None) -> dict:
    """为工具结果添加调试和执行信息
    
    Args:
        result: 原始结果字典
        cmd: 执行的命令
        params_file: 参数文件路径
        stdout: 标准输出
        stderr: 标准错误
        execution_time: 执行时间
        additional_info: 额外信息
    
    Returns:
        增强后的结果字典
    """
    try:
        # 构建调试信息
        debug_info = {}
        
        if cmd:
            debug_info["command"] = cmd
        if params_file:
            debug_info["params_file"] = params_file
        if execution_time > 0:
            debug_info["execution_time_seconds"] = round(execution_time, 2)
        if stdout:
            debug_info["stdout_preview"] = stdout[:500] + "..." if len(stdout) > 500 else stdout
        if stderr:
            debug_info["stderr_preview"] = stderr[:500] + "..." if len(stderr) > 500 else stderr
        if additional_info:
            debug_info.update(additional_info)
        
        # 添加时间戳
        debug_info["debug_timestamp"] = datetime.now().isoformat()
        
        # 将调试信息添加到结果中
        if debug_info:
            result["debug_info"] = debug_info
        
        # 记录调试日志
        logger.debug(f"工具调试信息已添加: {len(debug_info)} 项")
    
    except Exception as e:
        logger.warning(f"添加调试信息失败: {e}")
    
    return result