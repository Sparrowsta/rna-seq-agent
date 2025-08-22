"""
检测装饰器模块 - 简化检测工具的重复代码
专注于纯数据收集，不包含决策逻辑
"""

import os
import json
import time
import subprocess
from pathlib import Path
from functools import wraps
from typing import Dict, List, Any, Callable


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


def with_fastq_scan(func: Callable) -> Callable:
    """FASTQ文件扫描数据注入装饰器"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        fastq_data = _scan_fastq_files()
        return func(fastq_data, *args, **kwargs)
    return wrapper


def with_genome_config(func: Callable) -> Callable:
    """基因组配置数据注入装饰器"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        genome_config = _load_genome_config()
        return func(genome_config, *args, **kwargs)
    return wrapper


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


# ======================= 内部辅助函数 =======================

def _scan_fastq_files() -> Dict[str, Any]:
    """纯FASTQ文件扫描 - 只收集原始文件信息，不做任何判断"""
    project_root = Path(".")
    fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    
    # 扫描所有FASTQ文件
    all_fastq_files = []
    for ext in fastq_extensions:
        all_fastq_files.extend(project_root.rglob(ext))
    
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
    
    # 只返回原始文件列表，不做任何分析
    return {
        "detection_status": "success",
        "total_files_found": len(file_list),
        "fastq_files": file_list,
        "scan_timestamp": time.time()
    }


def _load_genome_config() -> Dict[str, Any]:
    """纯基因组配置加载 - 只检查文件存在性，不做状态判断"""
    genomes_file = Path("config/genomes.json")
    
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


def _extract_sample_name(filename: str) -> tuple:
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