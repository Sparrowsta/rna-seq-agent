"""
RNA-seq智能分析助手工具模块
提供数据收集功能，返回结构化数据供LLM智能处理和展示

重构原则:
- 工具专注纯数据收集，不做格式化
- 使用官方 @tool 装饰器
- 移除双模式逻辑，交由LLM决定展示方式
- 简化代码结构，提高维护性
"""

import json
import re
import time
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置系统
from .config import get_tools_config


# ==================== 工具函数 (使用 @tool 装饰器) ====================

@tool
def scan_fastq_files() -> Dict[str, Any]:
    """扫描项目中的FASTQ测序文件，返回文件列表、样本信息和基本统计数据"""
    config = get_tools_config()
    project_root = config.project_root
    fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    
    # 定义要排除的目录（中间文件和缓存目录）
    exclude_directories = {
        "work", "tmp", "temp", "results", "output", 
        ".nextflow", "logs", "cache", "__pycache__"
    }
    
    # 扫描所有FASTQ文件
    all_fastq_files = []
    for ext in fastq_extensions:
        for file_path in project_root.rglob(ext):
            if any(excluded in file_path.parts for excluded in exclude_directories):
                continue
            all_fastq_files.append(file_path)
    
    # 收集文件信息
    file_list = []
    for file_path in all_fastq_files:
        if not file_path.exists() or not file_path.is_file():
            continue
            
        try:
            file_size = file_path.stat().st_size
            
            # 提取样本信息
            filename = file_path.name
            if "_1." in filename or "_R1" in filename:
                sample_name = filename.split("_1.")[0].split("_R1")[0]
                read_type = "R1"
            elif "_2." in filename or "_R2" in filename:
                sample_name = filename.split("_2.")[0].split("_R2")[0]
                read_type = "R2"
            else:
                sample_name = filename.split(".")[0]
                read_type = "single"
            
            file_info = {
                "filename": filename,
                "full_path": str(file_path),
                "directory": str(file_path.parent),
                "size_bytes": file_size,
                "size_mb": round(file_size / 1024 / 1024, 2),
                "extension": file_path.suffix,
                "sample_name": sample_name,
                "read_type": read_type
            }
            file_list.append(file_info)
        except Exception:
            continue
    
    # 按文件名排序
    file_list.sort(key=lambda x: x["filename"])
    
    # 样本统计
    samples = {}
    for file_info in file_list:
        sample_name = file_info["sample_name"]
        read_type = file_info["read_type"]
        
        if sample_name not in samples:
            samples[sample_name] = {"R1": None, "R2": None, "single": None}
        
        samples[sample_name][read_type] = file_info
    
    # 确定测序类型
    paired_count = sum(1 for s in samples.values() if s["R1"] and s["R2"])
    single_count = sum(1 for s in samples.values() if s["single"])
    
    if paired_count > 0 and single_count == 0:
        sequencing_type = "paired_end"
    elif single_count > 0 and paired_count == 0:
        sequencing_type = "single_end"
    else:
        sequencing_type = "mixed"
    
    return {
        "detection_status": "success",
        "total_files": len(file_list),
        "total_samples": len(samples),
        "sequencing_type": sequencing_type,
        "paired_samples": paired_count,
        "single_samples": single_count,
        "total_size_mb": sum(f["size_mb"] for f in file_list),
        "files": file_list,
        "samples": samples,
        "scan_timestamp": time.time()
    }


@tool
def scan_system_resources() -> Dict[str, Any]:
    """检测系统硬件资源，返回CPU、内存、磁盘和负载信息"""
    try:
        import psutil
        
        # CPU信息
        cpu_count = psutil.cpu_count(logical=False) or 1
        logical_count = psutil.cpu_count(logical=True) or 1
        cpu_freq = psutil.cpu_freq()
        
        # 内存信息
        memory = psutil.virtual_memory()
        memory_gb = memory.total / 1024**3
        available_gb = memory.available / 1024**3
        used_percent = memory.percent
        
        # 磁盘信息
        disk = psutil.disk_usage('.')
        disk_total_gb = disk.total / 1024**3
        disk_free_gb = disk.free / 1024**3
        disk_used_percent = (disk.used / disk.total) * 100
        
        # 系统负载
        load_info = {}
        try:
            load_avg = psutil.getloadavg()
            load_ratio = load_avg[0] / max(logical_count, 1)
            load_info = {
                "load_1min": round(load_avg[0], 2),
                "load_5min": round(load_avg[1], 2),
                "load_15min": round(load_avg[2], 2),
                "load_ratio": round(load_ratio, 2)
            }
        except (AttributeError, OSError):
            load_info = {"error": "platform_not_supported"}
        
        return {
            "detection_status": "success",
            "cpu": {
                "physical_cores": cpu_count,
                "logical_cores": logical_count,
                "frequency_mhz": cpu_freq.current if cpu_freq else None
            },
            "memory": {
                "total_gb": round(memory_gb, 1),
                "available_gb": round(available_gb, 1),
                "used_percent": round(used_percent, 1),
                "total_bytes": memory.total,
                "available_bytes": memory.available
            },
            "disk": {
                "total_gb": round(disk_total_gb, 1),
                "free_gb": round(disk_free_gb, 1),
                "used_percent": round(disk_used_percent, 1),
                "total_bytes": disk.total,
                "free_bytes": disk.free
            },
            "load": load_info,
            "timestamp": time.time()
        }
        
    except ImportError:
        return {
            "detection_status": "missing_dependency",
            "error": "psutil not installed",
            "install_command": "uv add psutil"
        }
    except Exception as e:
        return {
            "detection_status": "error",
            "error": str(e)
        }


@tool
def scan_genome_files(genome_id: Optional[str] = None) -> Dict[str, Any]:
    """扫描可用的参考基因组配置，返回基因组列表和文件状态
    
    Args:
        genome_id: 可选的特定基因组ID，用于重点关注
    """
    config = get_tools_config()
    genomes_file = config.genomes_config_path
    
    if not genomes_file.exists():
        result = {"detection_status": "no_config_file"}
    else:
        try:
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
            
            if not genomes_data:
                result = {"detection_status": "empty_config"}
            else:
                # 检查每个基因组的文件状态
                genome_status = {}
                for genome_id_key, info in genomes_data.items():
                    fasta_path = info.get('fasta_path', '')
                    gtf_path = info.get('gtf_path', '')
                    
                    # 检查文件存在性
                    fasta_exists = bool(fasta_path and Path(fasta_path).exists())
                    gtf_exists = bool(gtf_path and Path(gtf_path).exists())
                    
                    # 检查索引状态
                    star_index_exists = False
                    hisat2_index_exists = False
                    
                    if fasta_exists:
                        # STAR索引检查
                        star_index_dir = config.get_star_index_dir(Path(fasta_path))
                        if star_index_dir.exists():
                            star_index_files = list(star_index_dir.iterdir())
                            star_index_exists = len(star_index_files) > 0
                        
                        # HISAT2索引检查
                        hisat2_index_dir = config.get_hisat2_index_dir(Path(fasta_path))
                        if hisat2_index_dir.exists():
                            ht2_files = list(hisat2_index_dir.glob("*.ht2"))
                            hisat2_index_exists = len(ht2_files) > 0
                    
                    genome_status[genome_id_key] = {
                        "species": info.get('species', ''),
                        "version": info.get('version', ''),
                        "fasta_url": info.get('fasta_url', ''),
                        "gtf_url": info.get('gtf_url', ''),
                        "fasta_path": fasta_path,
                        "gtf_path": gtf_path,
                        "fasta_exists": fasta_exists,
                        "gtf_exists": gtf_exists,
                        "complete": fasta_exists and gtf_exists,
                        "star_index_exists": star_index_exists,
                        "hisat2_index_exists": hisat2_index_exists
                    }
                
                result = {
                    "detection_status": "success",
                    "total_genomes": len(genomes_data),
                    "available_genomes": len([g for g in genome_status.values() if g["complete"]]),
                    "genomes": genome_status,
                    "config_path": str(genomes_file)
                }
        except Exception as e:
            result = {
                "detection_status": "error",
                "error": str(e)
            }
    
    # 如果指定了特定基因组，添加标记
    if genome_id:
        result["requested_genome"] = genome_id
    
    return result


@tool
def get_project_overview() -> Dict[str, Any]:
    """获取项目整体状态概览，包括数据、基因组、系统资源和分析历史"""
    return {
        "fastq_data": scan_fastq_files(),
        "genome_status": scan_genome_files(),
        "system_resources": scan_system_resources(),
        "analysis_history": list_analysis_history(),
        "overview_timestamp": time.time()
    }


@tool
def list_analysis_history() -> Dict[str, Any]:
    """获取历史分析记录，返回分析时间、配置和结果信息"""
    config = get_tools_config()
    reports_dir = config.reports_dir
    
    if not reports_dir.exists():
        return {
            "detection_status": "no_history",
            "total_analyses": 0,
            "analyses": []
        }
    
    # 扫描时间戳格式的归档文件夹
    analyses = []
    for item in reports_dir.iterdir():
        if item.is_dir() and not item.name.startswith('.') and item.name != "latest":
            # 检查是否是时间戳格式 (YYYYMMDD_HHMMSS)
            if len(item.name) == 15 and item.name[8] == '_':
                try:
                    timestamp = datetime.strptime(item.name, "%Y%m%d_%H%M%S")
                    
                    # 检查归档内容
                    analysis_info = {
                        "timestamp_str": item.name,
                        "timestamp": timestamp.timestamp(),
                        "path": str(item),
                        "files": {}
                    }
                    
                    # 检查各类文件
                    for file_name in ["analysis_report.json", "analysis_summary.md", 
                                    "runtime_config.json", "execution_log.txt"]:
                        file_path = item / file_name
                        if file_path.exists():
                            analysis_info["files"][file_name] = {
                                "exists": True,
                                "size_bytes": file_path.stat().st_size
                            }
                        else:
                            analysis_info["files"][file_name] = {"exists": False}
                    
                    # 尝试读取配置信息
                    runtime_config = item / "runtime_config.json"
                    if runtime_config.exists():
                        try:
                            with open(runtime_config, 'r', encoding='utf-8') as f:
                                config_data = json.load(f)
                            analysis_info["config"] = config_data.get("nextflow_params", {})
                        except Exception:
                            analysis_info["config"] = {}
                    
                    # 计算总大小
                    total_size = sum(f.stat().st_size for f in item.rglob("*") if f.is_file())
                    analysis_info["total_size_bytes"] = total_size
                    
                    analyses.append(analysis_info)
                except ValueError:
                    continue
    
    # 按时间排序（最新在前）
    analyses.sort(key=lambda x: x["timestamp"], reverse=True)
    
    return {
        "detection_status": "success",
        "total_analyses": len(analyses),
        "latest_analysis": analyses[0] if analyses else None,
        "analyses": analyses
    }


@tool
def check_tool_availability(tool_name: str) -> Dict[str, Any]:
    """检测生物信息学工具的可用性
    
    Args:
        tool_name: 工具名称 (fastp, star, hisat2, featurecounts)
    """
    tool_configs = {
        "fastp": ("qc_env", ["fastp", "--version"]),
        "star": ("align_env", ["STAR", "--version"]), 
        "hisat2": ("align_env", ["hisat2", "--version"]),
        "featurecounts": ("quant_env", ["featureCounts", "-v"])
    }
    
    if tool_name.lower() not in tool_configs:
        return {
            "tool_name": tool_name,
            "error": f"未知工具: {tool_name}",
            "available_tools": list(tool_configs.keys()),
            "available": False
        }
    
    env_name, cmd = tool_configs[tool_name.lower()]
    
    # 执行工具检测
    detection_data = {
        "tool_name": tool_name,
        "environment": env_name,
        "command": cmd,
        "timestamp": time.time()
    }
    
    try:
        # 执行版本检测命令
        full_cmd = ['micromamba', 'run', '-n', env_name] + cmd
        result = subprocess.run(full_cmd, capture_output=True, text=True, timeout=15)
        
        detection_data.update({
            "command_executed": True,
            "return_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "available": result.returncode == 0
        })
            
    except subprocess.TimeoutExpired:
        detection_data.update({
            "command_executed": False,
            "error": "timeout",
            "timeout_seconds": 15,
            "available": False
        })
    except FileNotFoundError:
        detection_data.update({
            "command_executed": False,
            "error": "micromamba_not_found",
            "available": False
        })
    except Exception as e:
        detection_data.update({
            "command_executed": False,
            "error": str(e),
            "available": False
        })
    
    return detection_data


@tool
def add_genome_config(user_input: str) -> Dict[str, Any]:
    """智能解析用户输入并添加基因组配置
    
    Args:
        user_input: 包含基因组信息和URL的用户输入
    """
    try:
        if not user_input.strip():
            return {
                "success": False,
                "error": "请提供基因组信息，包含FASTA和GTF文件的URL"
            }
        
        # 使用LLM解析用户输入
        from .core import get_shared_llm
        
        llm = get_shared_llm()
        
        system_prompt = """你是基因组信息解析专家。用户会提供包含FASTA和GTF文件URL的文本，你需要提取信息并返回JSON格式：

{
  "genome_id": "基因组标识符(如hg38,mm39,ce11)",
  "species": "物种名称(如human,mouse,caenorhabditis_elegans)", 
  "version": "版本号(通常与genome_id相同)",
  "fasta_url": "FASTA文件完整URL",
  "gtf_url": "GTF文件完整URL"
}

要求：
1. 从URL中智能识别基因组版本和物种
2. 常见映射：hg->human, mm->mouse, ce->caenorhabditis_elegans, dm->drosophila, rn->rat
3. 如果无法确定，根据生物信息学常识合理推断
4. 只返回有效的JSON，不要其他解释"""

        # 解析用户输入
        response = llm.invoke(f"System: {system_prompt}\n\nHuman: 请解析：{user_input}")
        parsed_content = str(response.content).strip()
        
        # 提取JSON内容
        try:
            genome_info = json.loads(parsed_content)
        except json.JSONDecodeError:
            json_match = re.search(r'\{.*\}', parsed_content, re.DOTALL)
            if not json_match:
                return {"success": False, "error": f"LLM解析失败：{parsed_content}"}
            genome_info = json.loads(json_match.group())
        
        # 验证必需字段
        required_fields = ['genome_id', 'species', 'version', 'fasta_url', 'gtf_url']
        missing_fields = [field for field in required_fields if not genome_info.get(field)]
        if missing_fields:
            return {"success": False, "error": f"缺少字段：{missing_fields}"}
        
        # 构建配置
        genome_id = genome_info['genome_id']
        species = genome_info['species']
        version = genome_info['version']
        
        fasta_path = f"genomes/{species}/{version}/{version}.fa"
        gtf_path = f"genomes/{species}/{version}/{version}.gtf"
        
        new_genome_config = {
            "species": species,
            "version": version,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path,
            "fasta_url": genome_info['fasta_url'],
            "gtf_url": genome_info['gtf_url']
        }
        
        # 读取现有配置并添加
        config = get_tools_config()
        genomes_file = config.genomes_config_path
        
        if genomes_file.exists():
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
        else:
            genomes_data = {}
        
        # 检查重复
        if genome_id in genomes_data:
            return {
                "success": False,
                "error": f"基因组 {genome_id} 已存在",
                "existing_config": genomes_data[genome_id]
            }
        
        # 保存配置
        genomes_data[genome_id] = new_genome_config
        config.path_manager.ensure_directory(config.settings.config_dir)
        
        with open(genomes_file, 'w', encoding='utf-8') as f:
            json.dump(genomes_data, f, indent=2, ensure_ascii=False)
        
        return {
            "success": True,
            "genome_id": genome_id,
            "config": new_genome_config,
            "message": f"成功添加基因组配置：{genome_id}"
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"添加基因组配置时出错：{str(e)}"
        }


@tool
def get_help() -> Dict[str, Any]:
    """获取系统帮助信息和功能说明"""
    return {
        "system_name": "RNA-seq智能分析助手",
        "current_mode": "Normal模式 (项目信息中心)",
        "core_tools": [
            "项目概览 - 一键查看项目完整状态和健康度",
            "历史分析记录 - 浏览已完成的分析和可复用配置",
            "查看FASTQ文件 - 详细扫描所有测序数据文件",
            "查看基因组信息 - 显示可用参考基因组状态",
            "添加基因组配置 - 智能解析并添加新的参考基因组"
        ],
        "next_steps": [
            "首次使用建议运行 '项目概览' 了解项目状态",
            "输入 '/plan' 进入计划模式进行深度分析规划"
        ],
        "mode_description": "Normal模式专注快速信息查看和项目概览，Plan模式负责深度检测和分析方案制定"
    }