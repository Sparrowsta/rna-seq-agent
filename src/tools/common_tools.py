"""
RNA-seq智能分析助手 - 通用检测和扫描工具

包含：
- scan_fastq_files: FASTQ文件扫描
- scan_system_resources: 系统资源检测  
- scan_genome_files: 基因组文件状态检查
"""

import json
import time
import psutil
from pathlib import Path
from typing import Dict, Any, Optional

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.common")


@tool
def scan_fastq_files() -> Dict[str, Any]:
    """扫描FASTQ文件，优先在数据目录(data/fastq)下查找，兼容容器挂载目录。

    返回：文件列表、样本信息和基本统计数据。
    """
    config = get_tools_config()
    fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

    # 定义要排除的目录（中间文件和缓存目录）
    exclude_directories = {
        "work", "tmp", "temp", "results", "cache"
    }

    # 选择搜索根目录：优先 data/fastq，其次 data，最后项目根目录
    search_roots = []
    try:
        if config.fastq_dir.exists():
            search_roots.append(config.fastq_dir)
        elif config.settings.data_dir.exists():
            search_roots.append(config.settings.data_dir)
        else:
            search_roots.append(config.project_root)
    except Exception:
        search_roots.append(config.project_root)

    # 扫描所有FASTQ文件
    all_fastq_files = []
    for root in search_roots:
        for ext in fastq_extensions:
            for file_path in root.rglob(ext):
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
    
    result = {
        "detection_status": "success",
        "search_roots": [str(p) for p in search_roots],
        "total_files": len(file_list),
        "total_samples": len(samples),
        "sequencing_type": sequencing_type,
        "paired_samples": paired_count,
        "single_samples": single_count,
        "total_size_mb": sum(f["size_mb"] for f in file_list),
        "files": file_list,
        "samples": samples,
    }
    # 日志：摘要与告警/调试
    try:
        if result["total_files"] == 0:
            logger.warning(f"FASTQ扫描：未找到任何文件，roots={result['search_roots']}")
        else:
            logger.info(
                f"FASTQ扫描：files={result['total_files']} samples={result['total_samples']} "
                f"paired={result['paired_samples']} single={result['single_samples']} type={result['sequencing_type']}"
            )
            # 预览前5个样本名
            sample_names = list(result["samples"].keys())[:5]
            logger.debug(f"样本预览(前5)：{sample_names}")
    except Exception:
        pass
    return result


@tool
def scan_system_resources() -> Dict[str, Any]:
    """检测系统硬件资源，返回CPU、内存、磁盘和负载信息"""
    try:
        # CPU信息
        cpu_count = psutil.cpu_count(logical=True)
        cpu_physical = psutil.cpu_count(logical=False)
        
        # 内存信息
        memory = psutil.virtual_memory()
        memory_gb = round(memory.total / (1024**3), 2)
        memory_available_gb = round(memory.available / (1024**3), 2)
        
        # 磁盘信息（工作目录）
        disk = psutil.disk_usage('/')
        disk_total_gb = round(disk.total / (1024**3), 2)
        disk_free_gb = round(disk.free / (1024**3), 2)
        
        # 系统负载（1分钟平均）
        try:
            load_avg = psutil.getloadavg()[0]
        except (AttributeError, OSError):
            load_avg = None
        
        load_info = {
            "1min_avg": load_avg,
            "cpu_usage_percent": psutil.cpu_percent(interval=1),
        }
        
        result = {
            "detection_status": "success",
            "cpu": {
                "logical_cores": cpu_count,
                "physical_cores": cpu_physical,
            },
            "memory": {
                "total_gb": memory_gb,
                "available_gb": memory_available_gb,
                "used_percent": memory.percent
            },
            "disk": {
                "total_gb": disk_total_gb,
                "free_gb": disk_free_gb,
                "used_percent": round((disk.used / disk.total) * 100, 1),
                "total_bytes": disk.total,
                "free_bytes": disk.free
            },
            "load": load_info,
            "timestamp": time.time()
        }
        # 日志：资源摘要与低资源告警
        try:
            logger.info(f"系统资源：CPU={cpu_count}核 内存={memory_gb}GB 磁盘={disk_free_gb}GB可用")
            if memory_available_gb < 4:
                logger.warning(f"内存资源告警：可用内存仅{memory_available_gb}GB")
            if disk_free_gb < 10:
                logger.warning(f"磁盘空间告警：可用空间仅{disk_free_gb}GB")
        except Exception:
            pass
        return result
    except Exception as e:
        return {
            "detection_status": "failed",
            "error": str(e),
            "timestamp": time.time()
        }


@tool
def scan_genome_files(genome_id: Optional[str] = None) -> Dict[str, Any]:
    """扫描可用的参考基因组配置，返回基因组列表和文件状态
    
    Args:
        genome_id: 可选的特定基因组ID，提供时只检查该基因组状态
    """
    config = get_tools_config()
    
    # 读取基因组配置（统一从 settings.genomes_config_path/data/genomes.json 处获取）
    genomes_config_path = config.genomes_config_path
    
    if not genomes_config_path.exists():
        return {
            "detection_status": "no_config",
            "available_genomes": 0,
            "genomes": {},
            "missing_config_file": str(genomes_config_path)
        }
    
    try:
        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
    except Exception as e:
        return {
            "detection_status": "config_error",
            "error": str(e),
            "config_path": str(genomes_config_path)
        }
    
    # 检查基因组文件状态
    genome_statuses = {}
    
    # 如果指定了genome_id，只检查该基因组
    genomes_to_check = {genome_id: genomes_config[genome_id]} if genome_id and genome_id in genomes_config else genomes_config
    
    for gid, genome_info in genomes_to_check.items():
        try:
            status = {
                "genome_id": gid,
                "species": genome_info.get("species", "unknown"),
                "version": genome_info.get("version", "unknown"),
                "files": {},
                "available": True,
                "missing_files": []
            }
            
            # 检查FASTA文件
            fasta_path = genome_info.get("fasta_path", "")
            if fasta_path:
                fasta_full_path = config.project_root / fasta_path
                status["files"]["fasta"] = {
                    "path": str(fasta_full_path),
                    "relative_path": fasta_path,
                    "exists": fasta_full_path.exists(),
                    "size_mb": round(fasta_full_path.stat().st_size / 1024 / 1024, 2) if fasta_full_path.exists() else 0
                }
                if not fasta_full_path.exists():
                    status["missing_files"].append("fasta")
                    status["available"] = False
            
            # 检查GTF文件
            gtf_path = genome_info.get("gtf_path", "")
            if gtf_path:
                gtf_full_path = config.project_root / gtf_path
                status["files"]["gtf"] = {
                    "path": str(gtf_full_path),
                    "relative_path": gtf_path,
                    "exists": gtf_full_path.exists(),
                    "size_mb": round(gtf_full_path.stat().st_size / 1024 / 1024, 2) if gtf_full_path.exists() else 0
                }
                if not gtf_full_path.exists():
                    status["missing_files"].append("gtf")
                    status["available"] = False
            
            # 检查STAR索引
            star_index_dir = genome_info.get("star_index_dir", "")
            if star_index_dir:
                star_index_full_path = config.project_root / star_index_dir
                star_genome_file = star_index_full_path / "Genome"
                status["files"]["star_index"] = {
                    "index_dir": str(star_index_full_path),
                    "relative_path": star_index_dir,
                    "exists": star_genome_file.exists(),
                    "complete": star_genome_file.exists() and (star_index_full_path / "SA").exists()
                }
                if not status["files"]["star_index"]["complete"]:
                    status["missing_files"].append("star_index")
            
            # 检查HISAT2索引
            hisat2_index_dir = genome_info.get("hisat2_index_dir", "")
            if hisat2_index_dir:
                hisat2_index_full_path = config.project_root / hisat2_index_dir
                # HISAT2索引通常是多个.ht2文件
                ht2_files = list(hisat2_index_full_path.glob("*.ht2"))
                status["files"]["hisat2_index"] = {
                    "index_dir": str(hisat2_index_full_path),
                    "relative_path": hisat2_index_dir,
                    "exists": len(ht2_files) > 0,
                    "index_files_count": len(ht2_files)
                }
                if len(ht2_files) == 0:
                    status["missing_files"].append("hisat2_index")
            
            genome_statuses[gid] = status
        except Exception as e:
            genome_statuses[gid] = {
                "genome_id": gid,
                "error": str(e),
                "available": False
            }
    
    # 统计
    available_count = sum(1 for status in genome_statuses.values() if status.get("available", False))
    
    result = {
        "detection_status": "success",
        "available_genomes": available_count,
        "total_configured": len(genome_statuses),
        "genomes": genome_statuses,
        "config_file": str(genomes_config_path)
    }
    # 日志：基因组概览
    try:
        if available_count == 0:
            logger.warning("基因组扫描：无可用基因组")
        else:
            available_ids = [gid for gid, status in genome_statuses.items() if status.get("available")]
            logger.info(f"基因组扫描：{available_count}/{len(genome_statuses)}可用 = {available_ids}")
    except Exception:
        pass
    return result
