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
from typing import Dict, Any, List

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
                "available_genomes": genome_info.get("summary", {}).get("available_count", 0),
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
        
        logger.info(f"项目概览: {fastq_info.get('total_samples', 0)}个样本, {genome_info.get('summary', {}).get('available_count', 0)}个基因组可用")
        return overview
        
    except Exception as error:
        logger.error(f"获取项目概览失败: {error}")
        return {
            "detection_status": "failed",
            "error": str(error),
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
                completed_stages = sum(1 for value in analysis_info["files"].values() if value)
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
        
    except Exception as error:
        logger.error(f"获取分析历史失败: {error}")
        return {
            "detection_status": "failed",
            "error": str(error),
            "total_analyses": 0,
            "analyses": []
        }


def write_params_file(step: str, params: dict, state: AgentState = None,
                     results_dir: str = None, metadata: dict = None) -> Path:
    """写入Nextflow参数文件

    Args:
        step: 执行步骤名称 (prepare/fastp/star/featurecounts等)
        params: 参数字典，直接写入文件供Nextflow使用
        state: 当前Agent状态 (可选，向后兼容)
        results_dir: 结果目录路径 (可选，优先使用)
        metadata: 额外的元数据 (可选，暂未使用)

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
        
        # 准备写入的数据 - 直接写入参数供Nextflow使用
        data_to_write = params
        
        # 写入文件
        with open(params_file, 'w', encoding='utf-8') as file_handle:
            json.dump(data_to_write, file_handle, indent=2, ensure_ascii=False)
        
        logger.info(f"参数文件已写入: {params_file}")
        return params_file
        
    except Exception as error:
        logger.error(f"写入参数文件失败: {error}")
        # 返回默认路径，避免中断流程（在异常中安全确定临时目录）
        try:
            tools_config = get_tools_config()
            temporary_directory = Path(tools_config.settings.temp_dir)
        except Exception:
            temporary_directory = Path("/tmp")
        try:
            temporary_directory.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass
        fallback_path = temporary_directory / f"{step}_params_error.json"
        return fallback_path


def extract_genome_paths(state: AgentState) -> Dict[str, str]:
    """从 AgentState 中动态提取简化的基因组路径信息

    Args:
        state: Agent状态，包含nextflow_config和query_results

    Returns:
        简化的基因组路径信息字典，包含各种路径的直接引用
        {
            "genome_id": "hg38",
            "species": "human", 
            "version": "hg38",
            "gtf_path": "/data/genomes/human/hg38/annotation.gtf",
            "star_index_path": "/data/genomes/human/hg38/star_index",
            "hisat2_index_path": "/data/genomes/human/hg38/hisat2_index"
        }
    """
    try:
        # 获取目标基因组ID
        genome_version = state.nextflow_config.get("genome_version") if state.nextflow_config else None
        if not genome_version:
            logger.warning("未在nextflow_config中找到genome_version")
            return {}

        # 从query_results中提取基因组信息
        query_results = state.query_results or {}
        genome_scan_result = query_results.get("verify_genome_setup", {})
        genomes_dict = genome_scan_result.get("genomes", {})

        # 获取目标基因组信息
        genome_info = genomes_dict.get(genome_version, {})

        if not genome_info:
            logger.warning(f"未找到基因组信息: {genome_version}, 可用基因组: {list(genomes_dict.keys())}")
            return {}

        # 提取文件路径信息
        files_info = genome_info.get("files", {})
        
        # 构建简化的路径字典
        paths = {
            "genome_id": genome_info.get("genome_id", genome_version),
            "species": genome_info.get("species", ""),
            "version": genome_info.get("version", genome_version),
            "gtf_path": files_info.get("gtf", {}).get("path", ""),
            "star_index_path": files_info.get("star_index", {}).get("path", ""),
            "hisat2_index_path": files_info.get("hisat2_index", {}).get("path", "")
        }

        # 验证关键路径是否存在
        missing_paths = []
        if not paths["gtf_path"]:
            missing_paths.append("gtf_path")
        if not paths["star_index_path"]:
            missing_paths.append("star_index_path")
        if not paths["hisat2_index_path"]:
            missing_paths.append("hisat2_index_path")

        if missing_paths:
            logger.warning(f"基因组 {genome_version} 缺少路径信息: {missing_paths}")

        logger.info(f"成功提取基因组路径信息: {genome_version}")
        return paths

    except Exception as e:
        logger.error(f"提取基因组路径信息失败: {e}")
        return {}


def normalize_resources(stage_name: str, resource_config: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """归一化资源配置为Nextflow可用的格式（仅支持工具同名键）

    Args:
        stage_name: 工具/阶段名称（fastp/star/hisat2/featurecounts）
        resource_config: 状态中的资源配置字典，键必须为工具名

    Returns:
        归一化后的资源配置字典 {cpus: int, memory: str}
    """
    try:
        resource_config_map = resource_config or {}
        stage_config = resource_config_map.get(stage_name, {})
        
        # 提取并验证cpus参数
        cpus = None
        if stage_config.get("cpus"):
            try:
                cpus = int(stage_config["cpus"])
            except (ValueError, TypeError):
                logger.warning(f"无效的cpus值: {stage_config.get('cpus')}")
        
        # 提取并验证memory参数
        memory = stage_config.get("memory")
        if memory and not isinstance(memory, str):
            memory = str(memory)
        
        # 构建归一化结果（只包含有效值）
        normalized = {}
        if cpus is not None:
            normalized["cpus"] = cpus
        if memory:
            normalized["memory"] = memory
            
        logger.debug(f"归一化资源配置 {stage_name}: {normalized}")
        return normalized
        
    except Exception as e:
        logger.error(f"归一化资源配置失败 {stage_name}: {e}")
        return {}

        

