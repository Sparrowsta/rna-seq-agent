"""
RNA-seq智能分析助手工具模块
提供FASTQ文件查询、基因组管理、用户意图收集等核心功能

重构说明:
- 合并decorators.py中的代码，统一管理所有工具函数
- 按功能模块重新组织，让辅助函数和主要函数组织在一起
- 保持双模式工具架构不变
"""

import os
import json
import glob
import re
import time
import subprocess
from pathlib import Path
from functools import wraps
from typing import Dict, List, Any, Union, Callable, Optional


# ==================== 装饰器和辅助系统 ====================

def pure_detection(error_prefix: str = "检测"):
    """纯检测装饰器 - 统一检测工具的返回格式和错误处理
    
    Args:
        error_prefix: 错误消息前缀
    
    Returns:
        标准化的检测结果: {"result": str, "query_results": dict, "config_updates": {}}
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
                
                # 如果函数返回的已经是标准格式
                if isinstance(result, dict) and "query_results" in result:
                    return result
                
                # 如果函数返回的是查询数据，自动包装
                if isinstance(result, dict):
                    return {
                        "result": f"✅ {error_prefix}完成",
                        "query_results": result,
                        "config_updates": {}  # 检测不生成配置
                    }
                
                return {
                    "result": str(result),
                    "query_results": {},
                    "config_updates": {}
                }
                
            except Exception as e:
                return {
                    "result": f"❌ {error_prefix}时出错: {str(e)}",
                    "query_results": {"status": "error", "error": str(e)},
                    "config_updates": {}
                }
        return wrapper
    return decorator


def tool_detection(tool_name: str, environment: str, version_cmd: List[str]):
    """纯工具检测装饰器 - 只收集版本信息，不做可用性判断
    
    Args:
        tool_name: 工具名称
        environment: micromamba环境名
        version_cmd: 版本检测命令列表
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            detection_data = {
                "tool_name": tool_name,
                "environment": environment,
                "command": version_cmd
            }
            
            try:
                # 执行版本检测命令
                cmd = ['micromamba', 'run', '-n', environment] + version_cmd
                conda_result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
                
                detection_data.update({
                    "command_executed": True,
                    "return_code": conda_result.returncode,
                    "stdout": conda_result.stdout,
                    "stderr": conda_result.stderr
                })
                    
            except subprocess.TimeoutExpired:
                detection_data.update({
                    "command_executed": False,
                    "error": "timeout",
                    "timeout_seconds": 15
                })
            except FileNotFoundError:
                detection_data.update({
                    "command_executed": False,
                    "error": "micromamba_not_found"
                })
            except Exception as e:
                detection_data.update({
                    "command_executed": False,
                    "error": str(e)
                })
            
            return {
                "result": f"🔧 {tool_name}检测完成",
                "query_results": detection_data,
                "config_updates": {}
            }
        return wrapper
    return decorator


def get_system_info() -> Dict[str, Any]:
    """获取系统硬件信息"""
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
            "load": load_info
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


# ==================== FASTQ文件处理模块 ====================

def _scan_fastq_files() -> Dict[str, Any]:
    """纯FASTQ文件扫描 - 排除中间文件目录"""
    project_root = Path(".")
    fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    
    # 定义要排除的目录（中间文件和缓存目录）
    exclude_directories = {
        "work", "tmp", "temp", "results", "output", 
        ".nextflow", "logs", "cache", "__pycache__"
    }
    
    # 扫描所有FASTQ文件，但排除特定目录
    all_fastq_files = []
    for ext in fastq_extensions:
        for file_path in project_root.rglob(ext):
            # 检查文件路径是否包含排除的目录
            if any(excluded in file_path.parts for excluded in exclude_directories):
                continue
            all_fastq_files.append(file_path)
    
    # 收集每个文件的原始信息
    file_list = []
    for file_path in all_fastq_files:
        # 检查文件是否真实存在且可访问
        if not file_path.exists() or not file_path.is_file():
            continue
            
        try:
            file_size = file_path.stat().st_size
            file_info = {
                "filename": file_path.name,
                "full_path": str(file_path),
                "directory": str(file_path.parent),
                "size_bytes": file_size,
                "size_mb": round(file_size / 1024 / 1024, 2),
                "extension": file_path.suffix
            }
            file_list.append(file_info)
        except Exception:
            # 跳过无法访问的文件
            continue
    
    # 按文件名排序
    file_list.sort(key=lambda x: x["filename"])
    
    return {
        "detection_status": "success",
        "total_files_found": len(file_list),
        "fastq_files": file_list,
        "scan_timestamp": time.time()
    }


def _extract_sample_names(file_list: list) -> list:
    """从文件列表提取样本名称列表"""
    sample_names = set()
    for file_info in file_list:
        filename = file_info.get("filename", "")
        sample_name, _ = _extract_sample_name_and_type(filename)
        sample_names.add(sample_name)
    return list(sample_names)


def _extract_sample_name_and_type(filename: str) -> tuple:
    """提取样本名称和读取类型"""
    if "_1." in filename or "_R1" in filename:
        sample_name = filename.split("_1.")[0].split("_R1")[0]
        read_type = "R1"
    elif "_2." in filename or "_R2" in filename:
        sample_name = filename.split("_2.")[0].split("_R2")[0]
        read_type = "R2"
    else:
        sample_name = filename.split(".")[0]
        read_type = "single"
    
    return sample_name, read_type


def _analyze_sample_pairing(raw_data: dict) -> dict:
    """分析样本配对关系（用于detailed模式）"""
    # 这里可以添加更复杂的配对分析逻辑
    # 目前返回原始数据，后续可以扩展
    return raw_data


def _format_fastq_report(raw_data: dict, depth: str = "basic") -> str:
    """格式化FASTQ数据为用户友好的报告"""
    file_count = raw_data.get("total_files_found", 0)
    file_list = raw_data.get("fastq_files", [])
    
    if file_count == 0:
        return "在项目目录中未找到任何FASTQ文件"
    
    # 计算基本统计
    total_size_mb = sum(f["size_mb"] for f in file_list)
    
    # 按目录分组
    samples_by_dir = {}
    for file_info in file_list:
        directory = file_info["directory"]
        filename = file_info["filename"]
        
        if directory not in samples_by_dir:
            samples_by_dir[directory] = {}
        
        # 提取样本名称和读取类型
        sample_name, read_type = _extract_sample_name_and_type(filename)
        
        if sample_name not in samples_by_dir[directory]:
            samples_by_dir[directory][sample_name] = {"R1": None, "R2": None, "single": None}
        
        samples_by_dir[directory][sample_name][read_type] = {
            "filename": filename,
            "size_mb": file_info["size_mb"],
            "full_path": file_info["full_path"]
        }
    
    # 统计样本数量和类型
    total_samples = sum(len(samples) for samples in samples_by_dir.values())
    paired_count = sum(
        1 for samples in samples_by_dir.values()
        for files in samples.values()
        if files['R1'] and files['R2']
    )
    single_count = sum(
        1 for samples in samples_by_dir.values()
        for files in samples.values()
        if files['single']
    )
    
    # 生成报告
    result = f"📊 **FASTQ数据概览**\n\n"
    result += f"📈 **统计信息:**\n"
    result += f"   - 样本数量: {total_samples} 个\n"
    result += f"   - 数据大小: {total_size_mb:.1f} MB\n"
    result += f"   - 分布目录: {len(samples_by_dir)} 个\n"
    
    if paired_count > 0 and single_count == 0:
        result += f"   - 测序类型: 双端测序 (Paired-end)\n"
        sequencing_type = "双端测序"
    elif single_count > 0 and paired_count == 0:
        result += f"   - 测序类型: 单端测序 (Single-end)\n"
        sequencing_type = "单端测序"
    else:
        result += f"   - 测序类型: 混合类型 ({paired_count}双端 + {single_count}单端)\n"
        sequencing_type = "混合类型"
    
    # 估算读数和分析建议
    estimated_reads = total_size_mb * 4  # 粗略估算
    result += f"   - 预估读数: ~{estimated_reads:.0f}M reads\n"
    
    # 智能分析建议
    result += f"\n💡 **分析建议:**\n"
    if total_samples >= 3:
        result += f"   - 样本数量充足，适合差异表达分析\n"
    elif total_samples == 2:
        result += f"   - 最小比较组，可进行基础差异分析\n"
    else:
        result += f"   - 单样本，适合表达谱分析或质控\n"
    
    if sequencing_type == "双端测序":
        result += f"   - 双端数据质量较高，推荐标准流程\n"
    
    if estimated_reads >= 100:
        result += f"   - 测序深度充足，支持深度分析\n"
    elif estimated_reads >= 40:
        result += f"   - 测序深度适中，满足基础分析\n"
    else:
        result += f"   - 测序深度较低，注意质量控制\n"
    
    if depth == "detailed":
        # 详细样本信息
        result += f"\n📁 **详细样本信息:**\n"
        for directory, samples in samples_by_dir.items():
            result += f"\n📂 {directory}:\n"
            for sample_name, files in samples.items():
                if files["R1"] and files["R2"]:
                    total_sample_size = files['R1']['size_mb'] + files['R2']['size_mb']
                    result += f"   ✅ {sample_name}: 双端配对 ({total_sample_size:.1f}MB)\n"
                    result += f"      └─ R1: {files['R1']['filename']} ({files['R1']['size_mb']}MB)\n"
                    result += f"      └─ R2: {files['R2']['filename']} ({files['R2']['size_mb']}MB)\n"
                elif files["single"]:
                    result += f"   📄 {sample_name}: 单端文件 ({files['single']['size_mb']}MB)\n"
                    result += f"      └─ {files['single']['filename']}\n"
    
    # 添加后续步骤建议
    result += f"\n🚀 **后续步骤:**\n"
    result += f"   • 运行 '项目概览' 查看完整项目状态\n"
    result += f"   • 使用 '/plan' 开始配置分析流程\n"
    
    return result.strip()


# ==================== 基因组文件处理模块 ====================

def _load_genome_config() -> Dict[str, Any]:
    """纯基因组配置加载 - 只检查文件存在性，不做状态判断"""
    genomes_file = Path("/config/genomes.json")
    
    if not genomes_file.exists():
        return {"detection_status": "no_config_file"}
    
    try:
        with open(genomes_file, 'r', encoding='utf-8') as f:
            genomes_data = json.load(f)
        
        if not genomes_data:
            return {"detection_status": "empty_config"}
            
        # 只检查文件存在性，不做任何分类判断
        genome_files = {}
        for genome_id, info in genomes_data.items():
            fasta_path = info.get('fasta_path', '')
            gtf_path = info.get('gtf_path', '')
            
            # 检查文件信息
            fasta_info = None
            gtf_info = None
            star_index_info = None
            
            if fasta_path:
                fasta_file = Path(fasta_path)
                if fasta_file.exists():
                    fasta_info = {
                        "exists": True,
                        "path": fasta_path,
                        "size_bytes": fasta_file.stat().st_size,
                        "size_mb": round(fasta_file.stat().st_size / 1024**2, 1)
                    }
                else:
                    fasta_info = {"exists": False, "path": fasta_path}
            
            if gtf_path:
                gtf_file = Path(gtf_path)
                if gtf_file.exists():
                    gtf_info = {
                        "exists": True,
                        "path": gtf_path,
                        "size_bytes": gtf_file.stat().st_size,
                        "size_mb": round(gtf_file.stat().st_size / 1024**2, 1)
                    }
                else:
                    gtf_info = {"exists": False, "path": gtf_path}
            
            # 检查STAR索引目录
            if fasta_path and Path(fasta_path).exists():
                star_index_dir = Path(fasta_path).parent / "star_index"
                if star_index_dir.exists():
                    index_files = list(star_index_dir.iterdir())
                    star_index_info = {
                        "exists": True,
                        "path": str(star_index_dir),
                        "file_count": len(index_files)
                    }
                else:
                    star_index_info = {"exists": False, "path": str(star_index_dir)}
            
            genome_files[genome_id] = {
                "species": info.get('species', ''),
                "version": info.get('version', ''),
                "fasta_url": info.get('fasta_url', ''),
                "gtf_url": info.get('gtf_url', ''),
                "fasta_file": fasta_info,
                "gtf_file": gtf_info,
                "star_index": star_index_info
            }
        
        return {
            "detection_status": "success",
            "total_genomes": len(genomes_data),
            "genome_files": genome_files,
            "config_path": str(genomes_file)
        }
        
    except Exception as e:
        return {
            "detection_status": "error",
            "error": str(e)
        }


def _format_genome_report(genome_data: dict) -> str:
    """格式化基因组数据为用户友好的报告"""
    total_genomes = genome_data.get("total_genomes", 0)
    genome_files = genome_data.get("genome_files", {})
    highlighted_genome = genome_data.get("highlighted_genome")
    
    if total_genomes == 0:
        return "未找到基因组配置文件或配置为空"
    
    result = f"可用基因组 ({total_genomes} 个):\n\n"
    
    for genome_id, info in genome_files.items():
        species = info.get('species', '未知物种')
        version = info.get('version', genome_id)
        
        # 检查本地文件状态
        fasta_info = info.get('fasta_file', {})
        gtf_info = info.get('gtf_file', {})
        star_index_info = info.get('star_index', {}) or {}
        
        fasta_status = "✅ 已下载" if fasta_info.get('exists') else "❌ 未下载"
        gtf_status = "✅ 已下载" if gtf_info.get('exists') else "❌ 未下载"
        
        # 高亮显示特定基因组
        if highlighted_genome == genome_id:
            result += f"🎯 **{genome_id} ({species})** [重点关注]\n"
        else:
            result += f"🧬 {genome_id} ({species})\n"
            
        result += f"   - 版本: {version}\n"
        result += f"   - FASTA: {fasta_status}\n"
        result += f"   - GTF: {gtf_status}\n"
        
        if fasta_info.get('exists'):
            size_mb = fasta_info.get('size_mb', 0)
            result += f"   - FASTA大小: {size_mb} MB\n"
        
        if star_index_info.get('exists'):
            file_count = star_index_info.get('file_count', 0)
            result += f"   - STAR索引: ✅ {file_count}个文件\n"
        else:
            result += f"   - STAR索引: ❌ 未构建\n"
        
        result += "\n"
    
    return result.strip()


def _format_system_report(system_info: dict) -> str:
    """格式化系统信息为用户友好的报告"""
    if system_info.get("detection_status") == "error":
        return f"❌ 系统信息检测错误: {system_info.get('error', '')}"
    
    result = "💻 **系统硬件资源检测**\n\n"
    
    # 显示系统信息
    cpu_info = system_info.get("cpu", {})
    memory_info = system_info.get("memory", {})
    disk_info = system_info.get("disk", {})
    load_info = system_info.get("load", {})
    
    result += "🔧 **CPU资源:**\n"
    result += f"   - 物理核心: {cpu_info.get('physical_cores', 0)} 个\n"
    result += f"   - 逻辑核心: {cpu_info.get('logical_cores', 0)} 个\n"
    if cpu_info.get('frequency_mhz'):
        result += f"   - 基础频率: {cpu_info['frequency_mhz']:.0f} MHz\n"
    
    result += "\n🧠 **内存资源:**\n"
    result += f"   - 总内存: {memory_info.get('total_gb', 0):.1f} GB\n"
    result += f"   - 可用内存: {memory_info.get('available_gb', 0):.1f} GB\n"
    result += f"   - 内存使用率: {memory_info.get('used_percent', 0):.1f}%\n"
    
    result += "\n💾 **磁盘空间:**\n"
    result += f"   - 总容量: {disk_info.get('total_gb', 0):.1f} GB\n"
    result += f"   - 可用空间: {disk_info.get('free_gb', 0):.1f} GB\n"
    result += f"   - 使用率: {disk_info.get('used_percent', 0):.1f}%\n"
    
    if "error" not in load_info:
        result += "\n📊 **系统负载:**\n"
        result += f"   - 1分钟平均负载: {load_info.get('load_1min', 0):.2f}\n"
        result += f"   - 5分钟平均负载: {load_info.get('load_5min', 0):.2f}\n"
        result += f"   - 15分钟平均负载: {load_info.get('load_15min', 0):.2f}\n"
        result += f"   - 负载比率: {load_info.get('load_ratio', 0):.2f} (相对于CPU核心)\n"
    
    # 添加资源评估建议
    result += "\n💡 **资源评估:**\n"
    total_gb = memory_info.get('total_gb', 0)
    cpu_cores = cpu_info.get('physical_cores', 0)
    
    if total_gb >= 16:
        result += "   - 内存充足，支持大型基因组分析\n"
    elif total_gb >= 8:
        result += "   - 内存适中，适合中小型数据集\n"
    else:
        result += "   - 内存偏低，建议使用测试数据集\n"
    
    if cpu_cores >= 8:
        result += "   - CPU核心充足，支持并行处理\n"
    elif cpu_cores >= 4:
        result += "   - CPU核心适中，满足基本分析需求\n"
    else:
        result += "   - CPU核心较少，处理时间可能较长\n"
    
    return result.strip()


# ==================== 双模式工具（支持Normal和Detect节点） ====================

def scan_fastq_files(mode: str = "normal", depth: str = "basic") -> Union[str, dict]:
    """统一FASTQ文件扫描工具 - 支持双模式输出
    
    Args:
        mode: "normal" 返回格式化字符串供用户阅读
              "detect" 返回结构化dict供系统处理
        depth: "basic" 快速扫描基本信息
               "detailed" 深度分析配对关系
    
    Returns:
        根据mode返回字符串或字典
    """
    try:
        # 使用共享的底层扫描函数
        raw_data = _scan_fastq_files()
        
        if raw_data.get("detection_status") != "success":
            error_msg = f"FASTQ文件扫描失败: {raw_data.get('error', '未知错误')}"
            if mode == "detect":
                return {
                    "result": error_msg,
                    "query_results": raw_data,
                    "config_updates": {}
                }
            else:
                return error_msg
        
        if mode == "detect":
            # Detect模式：返回结构化数据供系统处理
            if depth == "detailed":
                # 为Detect节点增强配对分析
                enhanced_data = _analyze_sample_pairing(raw_data)
                simplified_data = {
                    "detection_status": enhanced_data["detection_status"],
                    "file_count": enhanced_data["total_files_found"],
                    "file_paths": [f["full_path"] for f in enhanced_data["fastq_files"]],
                    "samples": _extract_sample_names(enhanced_data["fastq_files"])
                }
            else:
                simplified_data = {
                    "detection_status": raw_data["detection_status"],
                    "file_count": raw_data["total_files_found"],
                    "file_paths": [f["full_path"] for f in raw_data["fastq_files"]],
                    "samples": _extract_sample_names(raw_data["fastq_files"])
                }
            
            return {
                "result": f"✅ 扫描完成，发现{raw_data['total_files_found']}个FASTQ文件",
                "query_results": simplified_data,
                "config_updates": {}
            }
        else:
            # Normal模式：返回格式化报告
            return _format_fastq_report(raw_data, depth)
            
    except Exception as e:
        error_msg = f"FASTQ文件扫描时出错: {str(e)}"
        if mode == "detect":
            return {
                "result": error_msg,
                "query_results": {"detection_status": "error", "error": str(e)},
                "config_updates": {}
            }
        else:
            return error_msg


def scan_genome_files(mode: str = "normal", genome_id: Optional[str] = None) -> Union[str, dict]:
    """统一基因组文件扫描工具 - 支持双模式输出
    
    Args:
        mode: "normal" 返回格式化字符串供用户阅读
              "detect" 返回结构化dict供系统处理
        genome_id: 特定基因组ID，None表示检查所有
    
    Returns:
        根据mode返回字符串或字典
    """
    try:
        # 使用共享的基因组配置加载函数
        genome_data = _load_genome_config()
        
        if genome_data.get("detection_status") not in ["success"]:
            error_msg = f"基因组配置检查失败: {genome_data.get('error', '配置文件问题')}"
            if mode == "detect":
                return {
                    "result": error_msg,
                    "query_results": genome_data,
                    "config_updates": {}
                }
            else:
                return error_msg
        
        # 如果指定了特定基因组，验证存在性并标记
        if genome_id and genome_data.get("genome_files"):
            if genome_id not in genome_data["genome_files"]:
                error_msg = f"未找到基因组配置: {genome_id}"
                if mode == "detect":
                    return {
                        "result": error_msg,
                        "query_results": {"detection_status": "error", "error": error_msg},
                        "config_updates": {}
                    }
                else:
                    return error_msg
            
            # 保持所有基因组信息，但标记特定基因组
            genome_data["highlighted_genome"] = genome_id
        
        if mode == "detect":
            # Detect模式：返回结构化数据
            # 简化genome_data以减少JSON复杂度
            genomes = genome_data.get("genome_files", {})
            available_genomes = []
            for name, info in genomes.items():
                is_complete = all([
                    info.get("fasta_file", {}).get("exists", False),
                    info.get("gtf_file", {}).get("exists", False)
                ])
                available_genomes.append({
                    "name": name,
                    "species": info.get("species"),
                    "complete": is_complete,
                    "has_star_index": (info.get("star_index") or {}).get("exists", False)
                })
            
            simplified_config = {
                "detection_status": genome_data["detection_status"],
                "total_genomes": len(genomes),
                "available_genomes": available_genomes
            }
            
            return {
                "result": f"✅ 检测到{len(genomes)}个基因组配置",
                "query_results": simplified_config,
                "config_updates": {}
            }
        else:
            # Normal模式：返回格式化报告
            return _format_genome_report(genome_data)
            
    except Exception as e:
        error_msg = f"基因组文件扫描时出错: {str(e)}"
        if mode == "detect":
            return {
                "result": error_msg,
                "query_results": {"detection_status": "error", "error": str(e)},
                "config_updates": {}
            }
        else:
            return error_msg


def scan_system_resources(mode: str = "normal") -> Union[str, dict]:
    """统一系统资源扫描工具 - 支持双模式输出
    
    Args:
        mode: "normal" 返回格式化字符串供用户阅读
              "detect" 返回结构化dict供系统处理
    
    Returns:
        根据mode返回字符串或字典
    """
    try:
        # 使用共享的系统信息获取函数
        system_info = get_system_info()
        
        if system_info.get("detection_status") == "missing_dependency":
            error_msg = "❌ 无法检测系统资源 (psutil未安装)\n   请安装: uv add psutil"
            if mode == "detect":
                return {
                    "result": error_msg,
                    "query_results": system_info,
                    "config_updates": {}
                }
            else:
                return error_msg
        
        if system_info.get("detection_status") == "error":
            error_msg = f"⚠️ 系统资源检测错误: {system_info.get('error', '')}"
            if mode == "detect":
                return {
                    "result": error_msg,
                    "query_results": system_info,
                    "config_updates": {}
                }
            else:
                return error_msg
        
        if mode == "detect":
            # Detect模式：返回结构化数据
            return {
                "result": "✅ 系统资源检测完成",
                "query_results": system_info,
                "config_updates": {}
            }
        else:
            # Normal模式：返回格式化报告
            return _format_system_report(system_info)
            
    except Exception as e:
        error_msg = f"系统资源扫描时出错: {str(e)}"
        if mode == "detect":
            return {
                "result": error_msg,
                "query_results": {"detection_status": "error", "error": str(e)},
                "config_updates": {}
            }
        else:
            return error_msg


# ==================== 生物信息学工具检测 ====================

@tool_detection("fastp", "qc_env", ["fastp", "--version"])
def check_fastp_availability() -> dict:
    """检测fastp工具可用性 - 使用装饰器"""
    pass


@tool_detection("STAR", "align_env", ["STAR", "--version"])
def check_star_availability() -> dict:
    """检测STAR工具可用性 - 使用装饰器"""
    pass


@tool_detection("featureCounts", "quant_env", ["featureCounts", "-v"])
def check_featurecounts_availability() -> dict:
    """检测featureCounts工具可用性 - 使用装饰器"""
    pass


# ==================== Normal模式专用工具 ====================

def get_project_overview(query: str = "") -> str:
    """项目全貌概览 - 整合所有关键信息的智能仪表板"""
    try:
        result = "🎯 **项目概览仪表板**\n\n"
        
        # 1. 数据状态总览
        result += "📊 **数据状态:**\n"
        
        # 扫描FASTQ文件
        project_root = Path(".")
        fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        all_fastq_files = []
        for ext in fastq_extensions:
            all_fastq_files.extend(project_root.rglob(ext))
        
        # 过滤原始文件
        excluded_dirs = ["work", "results", "tmp"]
        processed_indicators = ["trimmed", "fastp", "cutadapt", "filtered", "processed", "qc"]
        raw_fastq_files = []
        
        for file_path in all_fastq_files:
            if not file_path.exists() or any(excluded_dir in file_path.parts for excluded_dir in excluded_dirs):
                continue
            filename_lower = file_path.name.lower()
            if not any(indicator in filename_lower for indicator in processed_indicators):
                raw_fastq_files.append(file_path)
        
        # 统计样本信息
        total_samples = 0
        total_size_mb = 0
        sequencing_type = "未检测到"
        
        if raw_fastq_files:
            # 简单样本计数和大小统计
            sample_names = set()
            for file_path in raw_fastq_files:
                filename = file_path.name
                total_size_mb += file_path.stat().st_size / 1024 / 1024
                
                if "_1." in filename or "_R1" in filename:
                    sample_name = filename.split("_1.")[0].split("_R1")[0]
                    sample_names.add(sample_name)
                    sequencing_type = "双端测序 (Paired-end)"
                elif "_2." in filename or "_R2" in filename:
                    sample_name = filename.split("_2.")[0].split("_R2")[0]
                    sample_names.add(sample_name)
                else:
                    sample_name = filename.split(".")[0]
                    sample_names.add(sample_name)
                    sequencing_type = "单端测序 (Single-end)"
            
            total_samples = len(sample_names)
        
        result += f"   - 样本数量: {total_samples} 个\n"
        result += f"   - 数据大小: {total_size_mb:.1f} MB\n"
        result += f"   - 测序类型: {sequencing_type}\n"
        
        # 2. 基因组状态
        result += "\n🧬 **基因组状态:**\n"
        genomes_file = Path("/config/genomes.json")
        ready_genomes = 0
        total_genomes = 0
        
        if genomes_file.exists():
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
            total_genomes = len(genomes_data)
            
            for genome_id, info in genomes_data.items():
                fasta_path = info.get('fasta_path', '')
                gtf_path = info.get('gtf_path', '')
                if fasta_path and gtf_path and Path(fasta_path).exists() and Path(gtf_path).exists():
                    ready_genomes += 1
        
        result += f"   - 可用基因组: {total_genomes} 个\n"
        result += f"   - 就绪基因组: {ready_genomes} 个\n"
        
        # 3. 历史分析
        result += "\n📈 **历史分析:**\n"
        results_dir = Path("data/results")
        analysis_count = 0
        latest_analysis = "无"
        
        if results_dir.exists():
            analysis_dirs = [d for d in results_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
            analysis_count = len(analysis_dirs)
            
            if analysis_dirs:
                # 找最新的分析
                latest_dir = max(analysis_dirs, key=lambda x: x.stat().st_mtime)
                latest_analysis = latest_dir.name
        
        result += f"   - 历史分析: {analysis_count} 次\n"
        result += f"   - 最新分析: {latest_analysis}\n"
        
        # 4. 项目健康度评估
        result += "\n💡 **项目评估:**\n"
        
        health_score = 0
        suggestions = []
        
        if total_samples > 0:
            health_score += 40
        else:
            suggestions.append("未检测到FASTQ数据文件")
        
        if ready_genomes > 0:
            health_score += 30
        else:
            suggestions.append("需要下载基因组参考文件")
        
        if total_samples >= 3:
            health_score += 20
            suggestions.append("样本数量充足，适合差异表达分析")
        elif total_samples >= 2:
            health_score += 10
            suggestions.append("样本数量较少，考虑增加重复")
        
        if sequencing_type == "双端测序 (Paired-end)":
            health_score += 10
            suggestions.append("双端测序数据，质量较高")
        
        result += f"   - 项目健康度: {health_score}/100\n"
        
        if health_score >= 80:
            result += "   - 状态: ✅ 项目就绪，可开始分析\n"
        elif health_score >= 60:
            result += "   - 状态: ⚠️ 基本就绪，建议检查配置\n"
        else:
            result += "   - 状态: ❌ 需要完善项目配置\n"
        
        # 5. 智能建议
        if suggestions:
            result += "\n🚀 **智能建议:**\n"
            for suggestion in suggestions[:3]:  # 最多显示3个建议
                result += f"   • {suggestion}\n"
        
        result += "\n💡 输入 '/plan' 开始配置分析流程"
        
        return result.strip()
        
    except Exception as e:
        return f"生成项目概览时出错: {str(e)}"


def list_analysis_history(query: str = "") -> str:
    """历史分析管理 - 浏览和管理已完成的分析"""
    try:
        result = "📈 **分析历史记录**\n\n"
        
        results_dir = Path("data/results")
        if not results_dir.exists():
            return "📭 暂无分析历史记录\n\n💡 完成首次分析后，历史记录将显示在这里"
        
        # 扫描结果目录
        analysis_dirs = []
        for item in results_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                analysis_dirs.append(item)
        
        if not analysis_dirs:
            return "📭 results目录存在但无分析记录\n\n💡 运行分析后结果将保存在data/results/目录"
        
        # 按修改时间排序（最新在前）
        analysis_dirs.sort(key=lambda x: x.stat().st_mtime, reverse=True)
        
        result += f"📂 **发现 {len(analysis_dirs)} 个分析记录:**\n\n"
        
        for i, analysis_dir in enumerate(analysis_dirs):
            if i >= 10:  # 只显示前10个最新的分析
                break
                
            dir_name = analysis_dir.name
            modification_time = analysis_dir.stat().st_mtime
            
            # 转换时间戳为可读格式
            mod_time_str = time.strftime("%Y-%m-%d %H:%M", time.localtime(modification_time))
            
            result += f"📁 **{dir_name}**\n"
            result += f"   - 分析时间: {mod_time_str}\n"
            
            # 分析目录内容
            subdirs = []
            files = []
            total_size = 0
            
            try:
                for item in analysis_dir.iterdir():
                    if item.is_dir():
                        subdirs.append(item.name)
                    else:
                        files.append(item.name)
                        try:
                            total_size += item.stat().st_size
                        except:
                            pass
                
                total_size_mb = total_size / 1024 / 1024
                result += f"   - 结果大小: {total_size_mb:.1f} MB\n"
                
                # 识别分析类型
                analysis_types = []
                if any("fastp" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("质控")
                if any("star" in subdir.lower() or "hisat" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("比对")
                if any("bam" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("比对结果")
                if any("counts" in f.lower() or "feature" in f.lower() for f in files):
                    analysis_types.append("定量")
                if any("summary" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("报告")
                
                if analysis_types:
                    result += f"   - 分析步骤: {', '.join(analysis_types)}\n"
                else:
                    result += "   - 分析步骤: 未知或部分完成\n"
                
                # 检查关键结果文件
                key_files = []
                for f in files:
                    if f.endswith(('.html', '.json', '.md')):
                        key_files.append(f)
                
                if key_files:
                    result += f"   - 关键文件: {', '.join(key_files[:3])}\n"
                    if len(key_files) > 3:
                        result += f"     (还有{len(key_files)-3}个文件...)\n"
                
                # 评估分析完整性
                if "summary" in subdirs or any("report" in f.lower() for f in files):
                    result += "   - 状态: ✅ 分析完整\n"
                elif len(subdirs) >= 2:  # 至少有2个处理步骤
                    result += "   - 状态: ⚠️ 分析部分完成\n"
                else:
                    result += "   - 状态: ❌ 分析可能未完成\n"
                    
            except Exception as e:
                result += f"   - 状态: ❌ 目录访问错误: {str(e)}\n"
            
            result += "\n"
        
        # 添加统计信息
        if len(analysis_dirs) > 10:
            result += f"⏭️ 只显示了最新的10个分析，共有{len(analysis_dirs)}个历史记录\n\n"
        
        # 分析成功配置提取
        result += "🔄 **可复用配置:**\n"
        successful_configs = []
        
        for analysis_dir in analysis_dirs[:5]:  # 检查最新5个分析
            # 查找配置文件
            config_files = []
            for item in analysis_dir.rglob("*config*"):
                if item.is_file() and item.suffix in ['.json', '.yaml', '.yml']:
                    config_files.append(item)
            
            if config_files:
                successful_configs.append({
                    'name': analysis_dir.name,
                    'time': time.strftime("%m-%d %H:%M", time.localtime(analysis_dir.stat().st_mtime)),
                    'configs': len(config_files)
                })
        
        if successful_configs:
            for config in successful_configs:
                result += f"   • {config['name']} ({config['time']}) - {config['configs']}个配置文件\n"
            result += "\n💡 这些配置可以在Plan模式中复用\n"
        else:
            result += "   • 暂无可识别的配置文件\n"
        
        result += "\n🚀 输入 '/plan' 开始新的分析流程"
        
        return result.strip()
        
    except Exception as e:
        return f"获取分析历史时出错: {str(e)}"


def add_genome_config(user_input: str = "") -> str:
    """使用LLM智能解析用户输入并添加基因组配置
    
    Args:
        user_input: 用户输入，包含"添加基因组"和URL信息
    
    Returns:
        添加结果的文本描述
    """
    try:
        print(f"🛠️ add_genome_config被调用，用户输入: {user_input}")
        
        if not user_input.strip():
            return "请提供基因组信息，包含FASTA和GTF文件的URL"
        
        # 使用LLM解析用户输入
        from src.core import get_shared_llm
        
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
4. 只返回有效的JSON，不要其他解释

示例：
输入："添加基因组 https://path/ce11/ce11.fa.gz https://path/ce11/ce11.gtf.gz"
输出：{"genome_id":"ce11","species":"caenorhabditis_elegans","version":"ce11","fasta_url":"https://path/ce11/ce11.fa.gz","gtf_url":"https://path/ce11/ce11.gtf.gz"}"""

        # 使用简单的invoke方法
        response = llm.invoke(f"System: {system_prompt}\n\nHuman: 请解析：{user_input}")
        parsed_content = response.content.strip()
        
        print(f"🤖 LLM解析结果: {parsed_content}")
        
        # 提取JSON内容
        try:
            # 尝试直接解析
            genome_info = json.loads(parsed_content)
        except json.JSONDecodeError:
            # 如果失败，尝试提取JSON部分
            json_match = re.search(r'\{.*\}', parsed_content, re.DOTALL)
            if not json_match:
                return f"LLM解析失败，无法提取有效的JSON格式。原始响应：{parsed_content}"
            
            try:
                genome_info = json.loads(json_match.group())
            except json.JSONDecodeError as e:
                return f"JSON解析错误：{str(e)}\n原始响应：{parsed_content}"
        
        # 验证必需字段
        required_fields = ['genome_id', 'species', 'version', 'fasta_url', 'gtf_url']
        missing_fields = [field for field in required_fields if not genome_info.get(field)]
        if missing_fields:
            return f"LLM解析结果缺少必需字段：{missing_fields}\n解析结果：{genome_info}"
        
        print(f"✅ 解析成功: {genome_info}")
        
        # 执行实际的基因组添加逻辑
        genome_id = genome_info['genome_id']
        species = genome_info['species']
        version = genome_info['version']
        fasta_url = genome_info['fasta_url']
        gtf_url = genome_info['gtf_url']
        
        # 验证URL格式
        if not (fasta_url.startswith('http://') or fasta_url.startswith('https://')):
            return f"FASTA URL格式无效：{fasta_url}"
        if not (gtf_url.startswith('http://') or gtf_url.startswith('https://')):
            return f"GTF URL格式无效：{gtf_url}"
        
        # 生成本地路径
        fasta_path = f"data/genomes/{species}/{version}/{version}.fa"
        gtf_path = f"data/genomes/{species}/{version}/{version}.gtf"
        
        # 构建新的基因组配置
        new_genome_config = {
            "species": species,
            "version": version,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path,
            "fasta_url": fasta_url,
            "gtf_url": gtf_url
        }
        
        # 读取现有配置
        genomes_file = Path("/config/genomes.json")
        if genomes_file.exists():
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
        else:
            genomes_data = {}
        
        # 检查是否已存在
        if genome_id in genomes_data:
            existing = genomes_data[genome_id]
            return f"""基因组 {genome_id} 已存在

现有配置：
- 物种：{existing.get('species', '未知')}
- 版本：{existing.get('version', '未知')}
- FASTA URL：{existing.get('fasta_url', '未设置')}
- GTF URL：{existing.get('gtf_url', '未设置')}

如需更新，请使用不同的genome_id或先删除现有配置"""
        
        # 添加新配置
        genomes_data[genome_id] = new_genome_config
        
        # 保存配置
        os.makedirs("/config", exist_ok=True)
        with open(genomes_file, 'w', encoding='utf-8') as f:
            json.dump(genomes_data, f, indent=2, ensure_ascii=False)
        
        print(f"✅ 基因组配置添加成功: {genome_id}")
        
        return f"""✅ 成功添加基因组配置：{genome_id}

📋 基因组信息：
- 基因组ID：{genome_id}
- 物种：{species}
- 版本：{version}

📁 本地路径配置：
- FASTA：{fasta_path}
- GTF：{gtf_path}

🌐 下载源：
- FASTA URL：{fasta_url}
- GTF URL：{gtf_url}

💡 提示：文件将在首次运行分析时自动下载到指定路径"""
        
    except Exception as e:
        print(f"❌ add_genome_config出错: {str(e)}")
        return f"智能解析并添加基因组时出错：{str(e)}"


def get_help(query: str = "") -> str:
    """获取系统帮助信息"""
    return """
🎯 RNA-seq智能分析助手 - Normal模式 (项目信息中心)

📊 **核心项目工具:**
• 项目概览 - 一键查看项目完整状态和健康度
• 历史分析记录 - 浏览已完成的分析和可复用配置

📋 **详细信息查询:**
• 查看FASTQ文件 - 详细扫描所有测序数据文件
• 查看基因组信息 - 显示可用参考基因组状态

🗄️ **基因组管理:**
• 添加基因组配置 - 智能解析并添加新的参考基因组

🚀 **开始分析:**
输入 "/plan" 进入计划模式，Plan模式将执行深度数据检测、系统资源评估并制定智能分析方案

💡 **使用建议:**
1. 首次使用建议运行 "项目概览" 了解项目状态
2. 如需详细信息，使用具体的查询工具
3. 项目了解完成后，使用 "/plan" 进入深度分析规划

🔄 **模式分工:**
- Normal模式: 快速信息查看和项目概览
- Plan模式: 深度检测、就绪评估和分析方案制定
""".strip()