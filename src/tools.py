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

# 导入配置模块
from .config import get_tools_config
from .logging_bootstrap import get_logger
from .config.settings import Settings

logger = get_logger("rna.tools")




def enhance_tool_result_with_debug(result: dict, cmd: str = "", 
                                  params_file: str = "", stdout: str = "", 
                                  stderr: str = "") -> dict:
    """在DEBUG模式下增强工具结果，添加调试信息
    
    Args:
        result: 原始工具结果字典
        cmd: 执行的命令
        params_file: 参数文件路径  
        stdout: 标准输出
        stderr: 标准错误输出
        
    Returns:
        增强后的结果字典
    """
    settings = Settings()
    
    if settings.debug_mode:
        debug_info = {
            "debug_mode": True,
            "command": cmd,
            "params_file": params_file,
            "stdout": stdout[:2000] if stdout else "",  # 截断长输出
            "stderr": stderr[:1000] if stderr else ""   # 截断错误信息
        }
        result["debug_info"] = debug_info
        logger.debug(f"工具调试信息已添加: {len(debug_info)} 项")
    
    return result


# ==================== 工具函数 (使用 @tool 装饰器) ====================

@tool
def scan_fastq_files() -> Dict[str, Any]:
    """扫描FASTQ文件，优先在数据目录(data/fastq)下查找，兼容容器挂载目录。

    返回：文件列表、样本信息和基本统计数据。
    """
    config = get_tools_config()
    fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

    # 定义要排除的目录（中间文件和缓存目录）
    exclude_directories = {
        "work", "tmp", "temp", "results", "output",
        ".nextflow", "logs", "cache", "__pycache__"
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
        "scan_timestamp": time.time()
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
        import psutil
        
        # CPU信息
        cpu_count = psutil.cpu_count(logical=False) or 1
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
            load_ratio = load_avg[0] / max(cpu_count, 1)
            load_info = {
                "load_1min": round(load_avg[0], 2),
                "load_5min": round(load_avg[1], 2),
                "load_15min": round(load_avg[2], 2),
                "load_ratio": round(load_ratio, 2)
            }
        except (AttributeError, OSError):
            load_info = {"error": "platform_not_supported"}
        
        result = {
            "detection_status": "success",
            "cpu": {
                "physical_cores": cpu_count,
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
        # 日志：资源摘要与低资源告警
        try:
            logger.info(
                f"资源：CPU={cpu_count}核, Mem={result['memory']['available_gb']:.1f}/{result['memory']['total_gb']:.1f}GB可用, "
                f"DiskFree={result['disk']['free_gb']:.1f}GB/{result['disk']['total_gb']:.1f}GB"
            )
            if available_gb < 8 or disk_free_gb < 20:
                logger.warning(
                    f"低资源告警：可用内存={available_gb:.1f}GB, 可用磁盘={disk_free_gb:.1f}GB"
                )
        except Exception:
            pass
        return result
        
    except ImportError:
        logger.error("psutil 未安装，无法检测系统资源")
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
        try:
            logger.warning(f"基因组配置文件缺失: {genomes_file}")
        except Exception:
            pass
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
                    
                    star_index_dir_str = ""
                    hisat2_index_dir_str = ""
                    if fasta_exists:
                        # STAR索引检查
                        star_index_dir = config.get_star_index_dir(Path(fasta_path))
                        star_index_dir_str = str(star_index_dir)
                        if star_index_dir.exists():
                            star_index_files = list(star_index_dir.iterdir())
                            star_index_exists = len(star_index_files) > 0
                        
                        # HISAT2索引检查
                        hisat2_index_dir = config.get_hisat2_index_dir(Path(fasta_path))
                        hisat2_index_dir_str = str(hisat2_index_dir)
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
                        "star_index_dir": star_index_dir_str,
                        "hisat2_index_exists": hisat2_index_exists,
                        "hisat2_index_dir": hisat2_index_dir_str
                    }
                
                result = {
                    "detection_status": "success",
                    "total_genomes": len(genomes_data),
                    "available_genomes": len([g for g in genome_status.values() if g["complete"]]),
                    "genomes": genome_status,
                    "config_path": str(genomes_file)
                }
                try:
                    logger.info(
                        f"基因组配置：total={result['total_genomes']} available={result['available_genomes']} path={genomes_file}"
                    )
                    # 对缺失文件的基因组发出 warn
                    missing = [gid for gid, g in genome_status.items() if not g["complete"]]
                    if missing:
                        logger.warning(f"基因组缺失文件：{missing[:5]}{'...' if len(missing)>5 else ''}")
                except Exception:
                    pass
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
    # 使用 BaseTool.invoke 以避免在工具内部相互调用产生弃用警告
    return {
        "fastq_data": scan_fastq_files.invoke({}),
        "genome_status": scan_genome_files.invoke({}),
        "system_resources": scan_system_resources.invoke({}),
        "analysis_history": list_analysis_history.invoke({}),
        "overview_timestamp": time.time()
    }


@tool
def list_analysis_history() -> Dict[str, Any]:
    """获取历史分析记录，返回分析时间、配置和结果信息"""
    config = get_tools_config()
    reports_dir = config.reports_dir
    
    if not reports_dir.exists():
        logger.warning(f"历史分析目录不存在: {reports_dir}")
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
                    for file_name in ["analysis_report.json", "analysis_summary.md", "execution_log.txt"]:
                        file_path = item / file_name
                        if file_path.exists():
                            analysis_info["files"][file_name] = {
                                "exists": True,
                                "size_bytes": file_path.stat().st_size
                            }
                        else:
                            analysis_info["files"][file_name] = {"exists": False}
                    
                    # 计算总大小
                    total_size = sum(f.stat().st_size for f in item.rglob("*") if f.is_file())
                    analysis_info["total_size_bytes"] = total_size
                    
                    analyses.append(analysis_info)
                except ValueError:
                    continue
    
    # 按时间排序（最新在前）
    analyses.sort(key=lambda x: x["timestamp"], reverse=True)
    
    result = {
        "detection_status": "success",
        "total_analyses": len(analyses),
        "latest_analysis": analyses[0] if analyses else None,
        "analyses": analyses
    }
    try:
        latest = result["latest_analysis"]["timestamp_str"] if result["latest_analysis"] else None
        logger.info(f"历史分析：total={result['total_analyses']} latest={latest}")
    except Exception:
        pass
    return result


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
            
        # 日志
        if result.returncode == 0:
            logger.info(f"工具可用：{tool_name} (env={env_name})")
        else:
            stderr_snip = (result.stderr or '')[:200]
            logger.warning(f"工具不可用：{tool_name} rc={result.returncode} stderr={stderr_snip}")
        logger.debug(f"检测命令：{' '.join(full_cmd)}")
    except subprocess.TimeoutExpired:
        detection_data.update({
            "command_executed": False,
            "error": "timeout",
            "timeout_seconds": 15,
            "available": False
        })
        logger.warning(f"工具检测超时：{tool_name} (env={env_name}) 15s")
    except FileNotFoundError:
        detection_data.update({
            "command_executed": False,
            "error": "micromamba_not_found",
            "available": False
        })
        logger.error("micromamba 未找到，无法检测工具可用性")
    except Exception as e:
        detection_data.update({
            "command_executed": False,
            "error": str(e),
            "available": False
        })
        logger.error(f"工具检测异常：{tool_name} error={e}")
    
    return detection_data


@tool
def add_genome_config(genome_info: Dict[str, Any]) -> Dict[str, Any]:
    """添加基因组配置到 genomes.json（不在工具内调用LLM）。

    期望输入（由 Normal 节点的 LLM 已解析好并传入）：
    {
      "genome_id": "如 hg38/mm39/...",
      "species": "如 human/mouse/...",
      "version": "通常与 genome_id 相同",
      "fasta_url": "FASTA 完整下载 URL",
      "gtf_url": "GTF 完整下载 URL",
      // 可选：当需要自定义时可传入以下两项；否则自动按规范路径生成
      "fasta_path": "genomes/<species>/<version>/<version>.fa",
      "gtf_path": "genomes/<species>/<version>/<version>.gtf"
    }
    """
    try:
        # 兼容：若传入字符串，尝试按 JSON 解析
        if isinstance(genome_info, str):
            try:
                genome_info = json.loads(genome_info)
            except Exception:
                return {"success": False, "error": "genome_info 需为对象或可解析为JSON的字符串"}

        if not genome_info or not isinstance(genome_info, dict):
            return {"success": False, "error": "参数 genome_info 不能为空，且需为对象"}

        # 验证必需字段
        required_fields = ["genome_id", "species", "version", "fasta_url", "gtf_url"]
        missing_fields = [f for f in required_fields if not genome_info.get(f)]
        if missing_fields:
            return {"success": False, "error": f"缺少字段：{missing_fields}"}

        genome_id = str(genome_info["genome_id"]).strip()
        species = str(genome_info["species"]).strip()
        version = str(genome_info["version"]).strip()
        if not genome_id or not species or not version:
            return {"success": False, "error": "genome_id/species/version 不能为空"}

        # 规范化本地存放路径（允许外部覆盖）
        fasta_path = genome_info.get("fasta_path") or f"genomes/{species}/{version}/{version}.fa"
        gtf_path = genome_info.get("gtf_path") or f"genomes/{species}/{version}/{version}.gtf"

        new_genome_config = {
            "species": species,
            "version": version,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path,
            "fasta_url": genome_info["fasta_url"],
            "gtf_url": genome_info["gtf_url"]
        }

        # 读取现有配置并添加
        config = get_tools_config()
        genomes_file = config.genomes_config_path

        # 如文件不存在，初始化为空对象（尽量自愈）
        if genomes_file.exists():
            with open(genomes_file, "r", encoding="utf-8") as f:
                try:
                    genomes_data = json.load(f) or {}
                except json.JSONDecodeError:
                    genomes_data = {}
        else:
            genomes_file.parent.mkdir(parents=True, exist_ok=True)
            genomes_data = {}

        # 检查重复
        if genome_id in genomes_data:
            logger.warning(f"重复的基因组ID：{genome_id}")
            return {
                "success": False,
                "error": f"基因组 {genome_id} 已存在",
                "existing_config": genomes_data[genome_id]
            }

        # 保存配置
        genomes_data[genome_id] = new_genome_config

        with open(genomes_file, "w", encoding="utf-8") as f:
            json.dump(genomes_data, f, indent=2, ensure_ascii=False)

        logger.info(f"已添加基因组配置：{genome_id} -> {genomes_file}")
        return {
            "success": True,
            "genome_id": genome_id,
            "config": new_genome_config,
            "config_path": str(genomes_file),
            "message": f"成功添加基因组配置：{genome_id}"
        }

    except Exception as e:
        logger.error(f"添加基因组配置失败：{e}")
        return {"success": False, "error": f"添加基因组配置时出错：{str(e)}"}





# ==================== FastP专用工具函数 ====================

@tool
def run_nextflow_fastp(fastp_params: Dict[str, Any], sample_info: Dict[str, Any]) -> Dict[str, Any]:
    """执行Nextflow FastP质量控制流程
    
    Args:
        fastp_params: FastP参数字典，例如 {"qualified_quality_phred": 25, "length_required": 50}
        sample_info: 样本信息，包含sample_groups等
    
    Returns:
        执行结果字典，包含状态、输出路径、执行日志等
    """
    try:
        config = get_tools_config()
        
        # 验证必需参数
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

        # 结果目录（运行根目录）：优先使用 sample_info 提供；否则按时间戳生成到 data/results/<timestamp>
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = sample_info.get("results_dir") or str(config.settings.data_dir / "results" / f"{timestamp}")

        # 工作目录：统一到 /data/work
        run_id = results_dir.split('/')[-1] if '/' in str(results_dir) else timestamp
        work_dir = Path(base_data_path) / "work" / f"fastp_{run_id}"
        work_dir.mkdir(parents=True, exist_ok=True)
        results_dir = Path(results_dir)
        # 确保结果目录存在，避免 publishDir 目标不存在造成的发布失败
        try:
            results_dir.mkdir(parents=True, exist_ok=True)
            (results_dir / "fastp").mkdir(parents=True, exist_ok=True)

        except Exception:
            pass
        logger.info(f"FastP启动：样本组={len(sample_info.get('sample_groups', {}))} 运行目录: base={base_data_path} work={work_dir} results={results_dir}")
        
        nextflow_params: Dict[str, Any] = {}
        for raw_key, value in (fastp_params or {}).items():
            if value is None:
                continue
            key = str(raw_key)
            if key.startswith("--"):
                key = key[2:]
            if key == "thread":
                key = "threads"
            nextflow_params[key] = value
        
        # 添加样本组信息
        nextflow_params["sample_groups"] = sample_info["sample_groups"]
        
        # 设置结果目录和数据路径
        nextflow_params["results_dir"] = str(results_dir)
        nextflow_params["data"] = base_data_path
        
        # 创建Nextflow参数文件 - 保存到fastp子目录中
        fastp_dir = results_dir / "fastp"
        fastp_dir.mkdir(parents=True, exist_ok=True)
        params_file = fastp_dir / "fastp_params.json"
        with open(params_file, 'w', encoding='utf-8') as f:
            json.dump(nextflow_params, f, indent=2, ensure_ascii=False)
        
        # 构建Nextflow命令（兼容Docker与本地路径）
        # 1) 优先使用src/nextflow/目录下的 fastp.nf（本地开发）
        # 2) Docker镜像中 fastp.nf 位于 /src/nextflow/（见 Dockerfile COPY）
        nf_candidates = [
            config.settings.project_root / "src" / "nextflow" / "fastp.nf",
            Path("/src/nextflow/fastp.nf")
        ]
        nextflow_script = None
        for cand in nf_candidates:
            if cand.exists():
                nextflow_script = cand
                break
        if nextflow_script is None:
            return {
                "success": False,
                "error": "未找到 fastp.nf 脚本，请检查容器内是否存在 /src/nextflow/fastp.nf 或本地src/nextflow目录",
                "searched": [str(p) for p in nf_candidates]
            }

        cmd = [
            "nextflow", "run",
            str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
            "--data", str(config.settings.data_dir)
        ]
        # 记录可读命令串（供调试与回传）
        cmd_str = " ".join(cmd)

        logger.info("执行Nextflow FastP流水线...")
        logger.info(f"参数文件: {params_file}")
        # 正确显示并使用本函数创建的 Nextflow 工作目录
        logger.info(f"工作目录: {work_dir}")
        logger.info(f"结果目录: {results_dir}")
        
        # 执行Nextflow流水线
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30分钟超时
            cwd=config.settings.project_root
        )
        
        
        execution_time = time.time() - start_time
        
        if result.returncode == 0:
            # 解析输出结果
            sample_groups = sample_info.get("sample_groups", [])
            sample_count = len(sample_groups)
            
            # 基于约定的发布目录结构，构造每个样本的输出文件路径
            per_sample_outputs = []
            fastp_root = results_dir / "fastp"
            for item in sample_groups:
                sid = item.get("sample_id") or item.get("id")
                if not sid:
                    continue
                sdir = fastp_root / sid
                out = {
                    "sample_id": sid,
                    "html": str(sdir / f"{sid}.fastp.html"),
                    "json": str(sdir / f"{sid}.fastp.json"),
                }
                # 判断单双端
                r2 = item.get("read2")
                if r2:
                    out.update({
                        "trimmed_r1": str(sdir / f"{sid}_1.trimmed.fastq.gz"),
                        "trimmed_r2": str(sdir / f"{sid}_2.trimmed.fastq.gz"),
                        "is_paired": True,
                    })
                else:
                    out.update({
                        "trimmed_single": str(sdir / f"{sid}.single.trimmed.fastq.gz"),
                        "is_paired": False,
                    })
                per_sample_outputs.append(out)
            
            payload = {
                "success": True,
                "message": f"FastP质控完成，处理了{sample_count}个样本",
                "execution_time": execution_time,
                "results_dir": str(results_dir),
                # 返回正确的 Nextflow 工作目录，便于用户排查
                "work_dir": str(work_dir),
                "params_file": str(params_file),
                "sample_count": sample_count,
                "per_sample_outputs": per_sample_outputs
            }
            # 仅在调试模式返回详细日志
            try:
                if get_tools_config().settings.debug_mode:
                    payload.update({
                        "stdout": result.stdout,
                        "stderr": result.stderr,
                        "nextflow_params": nextflow_params
                    })
            except Exception:
                pass
            
            # 使用增强函数添加调试信息
            try:
                payload = enhance_tool_result_with_debug(
                    payload, cmd_str, str(params_file), result.stdout, result.stderr
                )
            except Exception:
                # 避免调试增强异常影响主流程
                logger.debug("FastP调试信息增强失败，忽略并返回主体结果")
            
            logger.info(f"FastP完成：samples={sample_count} results={results_dir}")
            return payload
        else:
            stderr_snip = (result.stderr or '')[:400]
            logger.warning(f"FastP失败：rc={result.returncode} results={results_dir} stderr={stderr_snip}")
            payload = {
                "success": False,
                "error": f"Nextflow执行失败 (返回码: {result.returncode})",
                "execution_time": execution_time,
            }
            try:
                if get_tools_config().settings.debug_mode:
                    payload.update({
                        "stdout": result.stdout,
                        "stderr": result.stderr,
                        "cmd": " ".join(cmd)
                    })
            except Exception:
                pass
            return payload
            
    except subprocess.TimeoutExpired:
        payload = {
            "success": False,
            "error": "Nextflow执行超时（30分钟）",
            "execution_time": time.time() - start_time
        }
        logger.warning("FastP执行超时：30分钟")
        return payload
    except Exception as e:
        payload = {
            "success": False,
            "error": f"执行FastP流水线时发生错误: {str(e)}",
            "execution_time": 0
        }
        logger.error(f"FastP异常：{e}")
        return payload


@tool  
def parse_fastp_results(results_directory: str) -> Dict[str, Any]:
    """解析FastP结果文件，提取客观质量指标供LLM分析
    
    Args:
        results_directory: FastP结果目录路径
    
    Returns:
        解析的质量指标字典，包含各样本的质量统计、过滤率等客观数据
        注意：此工具仅提供客观数据分析，不生成优化建议。优化建议由LLM基于这些数据智能生成。
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
                with open(json_file, 'r') as f:
                    fastp_data = json.load(f)
                
                # 提取样本ID（从文件名）
                sample_id = json_file.stem.replace('.fastp', '')
                
                # 提取质量指标
                summary = fastp_data.get("summary", {})
                before_filtering = summary.get("before_filtering", {})
                after_filtering = summary.get("after_filtering", {})
                
                # 读取数和碱基数
                reads_before = before_filtering.get("total_reads", 0)
                reads_after = after_filtering.get("total_reads", 0)
                bases_before = before_filtering.get("total_bases", 0)
                bases_after = after_filtering.get("total_bases", 0)
                
                # 质量指标
                q20_before = before_filtering.get("q20_rate", 0)
                q20_after = after_filtering.get("q20_rate", 0)
                q30_before = before_filtering.get("q30_rate", 0)
                q30_after = after_filtering.get("q30_rate", 0)
                
                # GC含量
                gc_before = before_filtering.get("gc_content", 0)
                gc_after = after_filtering.get("gc_content", 0)
                
                # 过滤率计算
                read_pass_rate = reads_after / reads_before if reads_before > 0 else 0
                base_pass_rate = bases_after / bases_before if bases_before > 0 else 0
                
                # 平均长度
                avg_length_before = bases_before / reads_before if reads_before > 0 else 0
                avg_length_after = bases_after / reads_after if reads_after > 0 else 0
                
                sample_metric = {
                    "sample_id": sample_id,
                    "json_file": str(json_file),
                    "reads_before": reads_before,
                    "reads_after": reads_after,
                    "bases_before": bases_before,
                    "bases_after": bases_after,
                    "read_pass_rate": round(read_pass_rate, 4),
                    "base_pass_rate": round(base_pass_rate, 4),
                    "q20_before": round(q20_before, 4),
                    "q20_after": round(q20_after, 4),
                    "q30_before": round(q30_before, 4),
                    "q30_after": round(q30_after, 4),
                    "gc_content_before": round(gc_before, 4),
                    "gc_content_after": round(gc_after, 4),
                    "avg_length_before": round(avg_length_before, 1),
                    "avg_length_after": round(avg_length_after, 1),
                    "quality_improvement": {
                        "q20_improvement": round(q20_after - q20_before, 4),
                        "q30_improvement": round(q30_after - q30_before, 4)
                    }
                }
                
                sample_metrics.append(sample_metric)
                
                # 累积总体统计
                overall_stats["total_reads_before"] += reads_before
                overall_stats["total_reads_after"] += reads_after
                overall_stats["total_bases_before"] += bases_before
                overall_stats["total_bases_after"] += bases_after
                overall_stats["q20_rates"].append(q20_after)
                overall_stats["q30_rates"].append(q30_after)
                overall_stats["gc_contents"].append(gc_after)
                
            except Exception as e:
                sample_metrics.append({
                    "sample_id": json_file.stem.replace('.fastp', ''),
                    "json_file": str(json_file),
                    "error": f"解析失败: {str(e)}"
                })
        
        # 计算总体指标
        total_samples = len([m for m in sample_metrics if "error" not in m])
        if total_samples > 0:
            overall_read_pass_rate = overall_stats["total_reads_after"] / overall_stats["total_reads_before"]
            overall_base_pass_rate = overall_stats["total_bases_after"] / overall_stats["total_bases_before"]
            avg_q20_rate = sum(overall_stats["q20_rates"]) / len(overall_stats["q20_rates"])
            avg_q30_rate = sum(overall_stats["q30_rates"]) / len(overall_stats["q30_rates"])
            avg_gc_content = sum(overall_stats["gc_contents"]) / len(overall_stats["gc_contents"])
        else:
            overall_read_pass_rate = 0
            overall_base_pass_rate = 0
            avg_q20_rate = 0
            avg_q30_rate = 0
            avg_gc_content = 0
        
        # 质量评估（仅提供客观指标，不生成优化建议）
        quality_assessment = {
            "overall_quality": "good" if avg_q30_rate > 0.85 else "moderate" if avg_q30_rate > 0.7 else "poor",
            "pass_rate_status": "good" if overall_read_pass_rate > 0.8 else "moderate" if overall_read_pass_rate > 0.6 else "poor"
        }
        
        result = {
            "success": True,
            "total_samples": total_samples,
            "results_directory": results_directory,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_reads_before": overall_stats["total_reads_before"],
                "total_reads_after": overall_stats["total_reads_after"],
                "total_bases_before": overall_stats["total_bases_before"],
                "total_bases_after": overall_stats["total_bases_after"],
                "overall_read_pass_rate": round(overall_read_pass_rate, 4),
                "overall_base_pass_rate": round(overall_base_pass_rate, 4),
                "average_q20_rate": round(avg_q20_rate, 4),
                "average_q30_rate": round(avg_q30_rate, 4),
                "average_gc_content": round(avg_gc_content, 4)
            },
            "quality_assessment": quality_assessment
        }
        try:
            logger.info(
                f"FastP结果: samples={result['total_samples']} pass_rate={result['overall_statistics']['overall_read_pass_rate']} avg_q30={result['overall_statistics']['average_q30_rate']}"
            )
            if sample_metrics:
                logger.debug(f"样本指标预览：{sample_metrics[0]}")
        except Exception:
            pass
        return result
        
    except Exception as e:
        logger.error(f"解析FastP结果失败：{e}")
        return {
            "success": False,
            "error": f"解析FastP结果失败: {str(e)}"
        }


# ==================== STAR工具函数 ====================

@tool
def download_genome_assets(genome_id: str) -> Dict[str, Any]:
    """下载指定基因组的FASTA和GTF文件
    
    Args:
        genome_id: 基因组标识（如"hg38", "mm39"），用于在genomes.json中查询下载信息
        force: 若目标文件已存在，是否强制覆盖下载（默认否）
    
    Returns:
        Dict: 下载结果信息
        {
            "success": bool,
            "fasta_path": str,      # 下载后的FASTA文件路径
            "gtf_path": str,        # 下载后的GTF文件路径
            "downloaded": List[str], # 实际下载的文件列表
            "skipped": List[str],   # 跳过的文件列表
            "errors": List[str]     # 错误信息列表
        }
    """
    try:
        tools_config = get_tools_config()
        
        # 读取基因组配置
        genomes_config_path = tools_config.genomes_config_path
        if not genomes_config_path.exists():
            return {
                "success": False,
                "error": f"基因组配置文件不存在: {genomes_config_path}"
            }
        
        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        
        if genome_id not in genomes_config:
            return {
                "success": False,
                "error": f"基因组 '{genome_id}' 在配置文件中不存在"
            }
        
        config = genomes_config[genome_id]
        fasta_url = config.get("fasta_url")
        gtf_url = config.get("gtf_url")
        fasta_relative = config.get("fasta_path")
        gtf_relative = config.get("gtf_path")
        
        if not all([fasta_url, gtf_url, fasta_relative, gtf_relative]):
            return {
                "success": False,
                "error": f"基因组 '{genome_id}' 配置不完整，缺少必要的URL或路径信息"
            }
        
        # 构建绝对路径
        project_root = tools_config.settings.project_root
        fasta_path = project_root / fasta_relative
        gtf_path = project_root / gtf_relative
        
        downloaded = []
        errors = []
        
        # 创建目录
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        gtf_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 准备下载任务
        download_info = [
            (fasta_url, fasta_path, "FASTA"),
            (gtf_url, gtf_path, "GTF")
        ]
        
        # 并行下载 - 同时启动两个curl进程
        processes = []
        temp_paths = []
        
        for url, target_path, file_type in download_info:
            # 确定临时文件路径
            if url.endswith('.gz'):
                temp_path = target_path.with_suffix('.gz.part')
            else:
                temp_path = target_path.with_suffix(target_path.suffix + '.part')
            temp_paths.append((temp_path, target_path, file_type))
            
            # 启动下载进程
            cmd = [
                "curl", "-L", "-fS", "--retry", "5", 
                "--retry-delay", "5", "--retry-connrefused",
                "-C", "-", "--progress-bar",
                "-o", str(temp_path), url
            ]
            process = subprocess.Popen(cmd)
            processes.append(process)
        
        # 等待所有下载完成
        download_results = [p.wait() for p in processes]
        
        # 处理下载结果和解压
        success_flags = [False, False]  # [FASTA, GTF]
        
        for i, (return_code, (temp_path, target_path, file_type)) in enumerate(zip(download_results, temp_paths)):
            if return_code != 0:
                errors.append(f"{file_type}下载失败，返回码: {return_code}")
                continue
            
            # 解压处理
            original_name = temp_path.name.rstrip('.part')
            if original_name.endswith('.gz'):
                # 使用gzip解压到目标文件
                try:
                    with open(target_path, 'w', encoding='utf-8') as output_file:
                        result = subprocess.run(
                            ["gzip", "-dc", str(temp_path)], 
                            stdout=output_file, 
                            stderr=subprocess.PIPE,
                            timeout=1800,
                            text=True
                        )
                        if result.returncode != 0:
                            errors.append(f"{file_type}解压失败")
                            continue
                except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
                    errors.append(f"{file_type}解压异常: {str(e)}")
                    continue
                
                # 检查解压结果
                if not target_path.exists() or target_path.stat().st_size == 0:
                    errors.append(f"{file_type}解压后文件为空")
                    continue
                
                # 删除临时文件
                temp_path.unlink()
            else:
                # 直接重命名
                temp_path.rename(target_path)
            
            downloaded.append(f"{file_type}: {target_path}")
            success_flags[i] = True  # 按索引标记成功
        
        success = all(success_flags)  # 所有文件都成功
        
        return {
            "success": success,
            "fasta_path": str(fasta_path),
            "gtf_path": str(gtf_path),
            "downloaded": downloaded,
            "errors": errors
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"下载基因组资源失败: {str(e)}"
        }


@tool
def build_star_index(
    genome_id: str,
    sjdb_overhang: Optional[int] = None,
    runThreadN: Optional[int] = None,
    force_rebuild: bool = False,
    results_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """构建STAR基因组索引（使用 build_index.nf）
    
    主要在RNA-seq分析流程中自动调用，当系统检测到缺少STAR索引时触发。
    参数文件统一保存到results/star目录，与其他STAR相关文件放在一起。

    Args:
        genome_id: 基因组标识，用于定位FASTA和GTF文件
        sjdb_overhang: 剪接位点overhang（默认用 DEFAULT_STAR_PARAMS 或 100）
        runThreadN: 线程数（默认用 DEFAULT_STAR_PARAMS 或 4）
        force_rebuild: 若索引目录已存在是否强制重建
        results_dir: 分析结果目录，参数文件保存到results/star/子目录

    Returns:
        Dict: 执行结果（包含 index_dir/stdout/stderr/skipped 等）
    """
    try:
        from .config.default_tool_params import DEFAULT_STAR_PARAMS
        tools_config = get_tools_config()

        # 读取基因组配置
        genomes_config_path = tools_config.genomes_config_path
        if not genomes_config_path.exists():
            return {"success": False, "error": f"基因组配置文件不存在: {genomes_config_path}"}

        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        if genome_id not in genomes_config:
            return {"success": False, "error": f"基因组 '{genome_id}' 在配置文件中不存在"}

        cfg = genomes_config[genome_id]
        fasta_rel = cfg.get("fasta_path")
        gtf_rel = cfg.get("gtf_path")
        if not fasta_rel or not gtf_rel:
            return {"success": False, "error": f"基因组 '{genome_id}' 缺少 fasta_path 或 gtf_path"}

        project_root = tools_config.settings.project_root
        fasta_path = project_root / fasta_rel
        gtf_path = project_root / gtf_rel
        
        # 检查文件是否存在，如果不存在则尝试下载
        files_missing = []
        files_invalid = []
        
        if not fasta_path.exists():
            files_missing.append(f"FASTA: {fasta_path}")
        else:
            # 检查FASTA文件有效性
            try:
                file_size = fasta_path.stat().st_size
                if file_size == 0:
                    files_invalid.append(f"FASTA文件为空: {fasta_path}")
                else:
                    # 检查文件格式（前几个字符）
                    with open(fasta_path, 'rb') as f:
                        first_bytes = f.read(10)
                    if not first_bytes.startswith(b'>'):
                        files_invalid.append(f"FASTA文件格式无效: {fasta_path} (首字符: {first_bytes[:5]})")
            except Exception as e:
                files_invalid.append(f"FASTA文件检查失败: {fasta_path} ({e})")
        
        if not gtf_path.exists():
            files_missing.append(f"GTF: {gtf_path}")
        else:
            # 检查GTF文件有效性
            try:
                file_size = gtf_path.stat().st_size
                if file_size == 0:
                    files_invalid.append(f"GTF文件为空: {gtf_path}")
            except Exception as e:
                files_invalid.append(f"GTF文件检查失败: {gtf_path} ({e})")
            
        if files_missing or files_invalid:
            logger.info(f"检测到问题文件，尝试重新下载: missing={files_missing}, invalid={files_invalid}")
            
            # 删除无效文件
            for invalid_file in files_invalid:
                try:
                    if "FASTA" in invalid_file and fasta_path.exists():
                        fasta_path.unlink()
                        logger.info(f"删除无效FASTA文件: {fasta_path}")
                    elif "GTF" in invalid_file and gtf_path.exists():
                        gtf_path.unlink()
                        logger.info(f"删除无效GTF文件: {gtf_path}")
                except Exception as e:
                    logger.warning(f"删除无效文件失败: {e}")
            
            download_result = download_genome_assets(genome_id)
            
            if not download_result.get("success", False):
                return {
                    "success": False,
                    "error": f"重新下载基因组文件失败: {download_result.get('error', '未知错误')}"
                }
            
            # 等待并再次检查文件
            import time
            time.sleep(2)  # 等待文件系统同步
            
            # 重新检查文件是否存在且有效
            if not fasta_path.exists():
                return {"success": False, "error": f"重新下载后FASTA文件仍不存在: {fasta_path}"}
            if not gtf_path.exists():
                return {"success": False, "error": f"重新下载后GTF文件仍不存在: {gtf_path}"}
            
            # 再次验证文件有效性
            try:
                fasta_size = fasta_path.stat().st_size
                if fasta_size == 0:
                    return {"success": False, "error": f"重新下载的FASTA文件仍为空: {fasta_path}"}
                
                with open(fasta_path, 'rb') as f:
                    first_bytes = f.read(10)
                if not first_bytes.startswith(b'>'):
                    return {"success": False, "error": f"重新下载的FASTA文件格式仍无效: {fasta_path}"}
                
                gtf_size = gtf_path.stat().st_size
                if gtf_size == 0:
                    return {"success": False, "error": f"重新下载的GTF文件仍为空: {gtf_path}"}
                    
            except Exception as e:
                return {"success": False, "error": f"验证重新下载的文件失败: {e}"}
            
            logger.info(f"文件重新下载并验证成功: FASTA={fasta_size}字节, GTF={gtf_size}字节")

        # 目标索引目录
        index_dir = tools_config.get_star_index_dir(fasta_path)

        # 存在即跳过
        if index_dir.exists() and not force_rebuild:
            key_files = ["SA", "SAindex", "Genome"]
            if all((index_dir / f).exists() for f in key_files):
                try:
                    logger.info(f"STAR索引已存在，跳过构建: {index_dir}")
                except Exception:
                    pass
                return {
                    "success": True,
                    "index_dir": str(index_dir),
                    "skipped": True,
                    "message": "STAR索引已存在，跳过构建",
                }

        # 构造 Nextflow 参数
        sjdb = (
            int(sjdb_overhang)
            if sjdb_overhang is not None
            else int(DEFAULT_STAR_PARAMS.get("sjdbOverhang") or 100)
        )
        threads = (
            int(runThreadN)
            if runThreadN is not None
            else int(DEFAULT_STAR_PARAMS.get("runThreadN", 4))
        )

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        temp_dir = Path(tools_config.settings.data_dir) / "tmp" / f"nextflow_star_index_{timestamp}"
        temp_dir.mkdir(parents=True, exist_ok=True)

        nf_params = {
            "genome_fasta": str(fasta_path),
            "genome_gtf": str(gtf_path),
            "star_index_dir": str(index_dir),
            "sjdb_overhang": sjdb,
            "runThreadN": threads,
            "limitGenomeGenerateRAM": 32000000000,
        }

        # 创建参数文件 - 统一保存到results/star目录（与star_params在一起）
        if results_dir:
            results_path = Path(results_dir)
            star_subdir = results_path / "star"
            star_subdir.mkdir(parents=True, exist_ok=True)
            params_file = star_subdir / "build_index_params.json"
        else:
            # 罕见的独立使用情况，仍然保存到基因组star目录
            star_dir = index_dir.parent / "star"
            star_dir.mkdir(parents=True, exist_ok=True)
            params_file = star_dir / f"build_index_params_{timestamp}.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 定位 build_index.nf
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "build_index.nf",
            Path("/src/nextflow/build_index.nf"),
        ]
        nextflow_script = None
        for cand in nf_candidates:
            if cand.exists():
                nextflow_script = cand
                break
        if nextflow_script is None:
            return {
                "success": False,
                "error": "未找到 build_index.nf 脚本，请检查 /src/nextflow/build_index.nf 路径",
                "searched": [str(p) for p in nf_candidates],
            }

        # 运行 Nextflow
        # 统一 Nextflow 工作目录到 /data/work
        work_root = tools_config.settings.data_dir / "work"
        work_dir = work_root / f"star_index_{timestamp}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # 调试日志：验证文件路径假设
        logger.info("=== STAR索引构建调试信息 ===")
        logger.info(f"FASTA文件路径: {fasta_path}")
        logger.info(f"FASTA文件存在: {fasta_path.exists()}")
        logger.info(f"GTF文件路径: {gtf_path}")
        logger.info(f"GTF文件存在: {gtf_path.exists()}")
        logger.info(f"目标索引目录: {index_dir}")
        logger.info(f"索引目录父目录存在: {index_dir.parent.exists()}")
        logger.info(f"Nextflow脚本路径: {nextflow_script}")
        logger.info(f"Nextflow脚本存在: {nextflow_script.exists()}")
        logger.info(f"参数文件路径: {params_file}")
        logger.info(f"工作目录: {work_dir}")
        
        # 调试日志：验证参数传递假设
        logger.info("=== Nextflow参数内容 ===")
        for key, value in nf_params.items():
            logger.info(f"{key}: {value}")
        
        # 验证参数文件内容
        try:
            with open(params_file, 'r', encoding='utf-8') as f:
                params_content = f.read()
            logger.info(f"参数文件内容: {params_content}")
        except Exception as e:
            logger.error(f"读取参数文件失败: {e}")
        
        cmd = [
            "nextflow",
            "run",
            str(nextflow_script),
            "-params-file",
            str(params_file),
            "-work-dir",
            str(work_dir),
        ]

        logger.info("构建STAR索引 (Nextflow)")
        logger.info(f"执行命令: {' '.join(cmd)}")
        logger.info(f"参数文件: {params_file}")
        logger.info(f"目标目录: {index_dir}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,
            cwd=tools_config.settings.project_root,
        )

        payload = {
            "success": result.returncode == 0,
            "index_dir": str(index_dir),
            "skipped": False,
            "params_file": str(params_file),
        }
        try:
            if get_tools_config().settings.debug_mode:
                payload.update({
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                    "cmd": " ".join(cmd)
                })
        except Exception:
            pass
        # 记录结果日志
        try:
            if payload["success"]:
                logger.info(f"STAR索引构建完成：{index_dir}")
            else:
                logger.warning(f"STAR索引构建失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        except Exception:
            pass
        return payload

    except subprocess.TimeoutExpired:
        try:
            logger.warning("STAR索引构建超时（120分钟）")
        except Exception:
            pass
        return {"success": False, "error": "STAR索引构建超时"}
    except Exception as e:
        try:
            logger.error(f"构建STAR索引异常：{e}")
        except Exception:
            pass
        return {"success": False, "error": f"构建STAR索引失败: {str(e)}"}


@tool
def run_nextflow_star(
    star_params: Dict[str, Any],
    fastp_results: Dict[str, Any],
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None
) -> Dict[str, Any]:
    """执行 STAR 比对（精简版）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - STAR 索引优先使用 genome_info.star_index_dir；否则由 genome_info.fasta_path 推导
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 star.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行STAR比对"}

        fastp_results_dir = fastp_results.get("results_dir")
        if not fastp_results_dir:
            return {"success": False, "error": "FastP结果缺少results_dir"}

        per_sample = fastp_results.get("per_sample_outputs") or []
        if not per_sample:
            return {"success": False, "error": "FastP结果缺少per_sample_outputs"}

        # 2) 运行根目录与工作目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = Path(fastp_results_dir)
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"star_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"STAR启动：samples={len(per_sample)} results={results_dir} work={work_dir}")

        # 3) 解析 STAR 索引目录
        def _resolve(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (tools_config.settings.project_root / pp)

        star_index_dir = ""
        if isinstance(genome_info, dict):
            star_index_dir = genome_info.get("star_index_dir") or genome_info.get("index_dir") or ""
            if not star_index_dir:
                fasta_path = genome_info.get("fasta_path") or genome_info.get("fasta")
                if fasta_path:
                    star_index_dir = str(tools_config.get_star_index_dir(_resolve(fasta_path)))

        if not star_index_dir:
            return {"success": False, "error": "缺少STAR索引目录（genome_info.star_index_dir 或 fasta_path 必须提供）"}

        star_index_path = _resolve(star_index_dir)
        if not star_index_path.exists():
            logger.warning(f"STAR索引不存在: {star_index_path}")
            return {"success": False, "error": f"STAR索引不存在: {star_index_path}"}

        # 4) 构造 sample_inputs（仅使用 FastP 返回的结构）
        sample_inputs: List[Dict[str, Any]] = []
        for i, fp in enumerate(per_sample):
            sid = fp.get("sample_id", f"sample_{i+1}")
            r1 = fp.get("trimmed_single") or fp.get("trimmed_r1")
            r2 = fp.get("trimmed_r2")
            if not r1:
                continue
            sample_inputs.append({
                "sample_id": sid,
                "is_paired": bool(r2),
                "read1": r1,
                **({"read2": r2} if r2 else {})
            })
        if not sample_inputs:
            return {"success": False, "error": "未从FastP结果构造到任何样本输入"}

        # 5) 组装 Nextflow 参数
        cleaned_params: Dict[str, Any] = {}
        for k, v in (star_params or {}).items():
            if v is None or k in {"star_cpus", "outFileNamePrefix"}:
                continue
            cleaned_params[k.lstrip('-')] = v

        nf_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            "star_index": str(star_index_path),
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **cleaned_params,
        }

        # 保存参数文件到star子目录
        star_dir = results_dir / "star"
        star_dir.mkdir(parents=True, exist_ok=True)
        params_file = star_dir / "star_params.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "star.nf",
            Path("/src/nextflow/star.nf"),
        ]
        nextflow_script = next((p for p in nf_candidates if p.exists()), None)
        if nextflow_script is None:
            return {"success": False, "error": "未找到 star.nf 脚本，请检查 /src/nextflow/star.nf 路径", "searched": [str(p) for p in nf_candidates]}

        logger.info(f"执行STAR比对 - 参数文件: {params_file}")
        logger.info(f"STAR索引: {nf_params['star_index']}")
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 star.nf publishDir 对齐）
        star_out = results_dir / "star"
        per_sample_outputs: List[Dict[str, Any]] = []
        for item in sample_inputs:
            sid = item["sample_id"]
            sdir = star_out / sid
            entry = {
                "sample_id": sid,
                "aligned_bam": str(sdir / f"{sid}.Aligned.sortedByCoord.out.bam"),
                "log_final": str(sdir / f"{sid}.Log.final.out"),
                "log_out": str(sdir / f"{sid}.Log.out"),
                "log_progress": str(sdir / f"{sid}.Log.progress.out"),
                "splice_junctions": str(sdir / f"{sid}.SJ.out.tab"),
            }
            qm = str(nf_params.get("quantMode", ""))
            if "TranscriptomeSAM" in qm:
                entry["transcriptome_bam"] = str(sdir / f"{sid}.Aligned.toTranscriptome.out.bam")
            if "GeneCounts" in qm:
                entry["gene_counts"] = str(sdir / f"{sid}.ReadsPerGene.out.tab")
            per_sample_outputs.append(entry)

        payload = {
            "success": result.returncode == 0,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": len(sample_inputs),
            "per_sample_outputs": per_sample_outputs,
        }
        if get_tools_config().settings.debug_mode:
            payload.update({"stdout": result.stdout, "stderr": result.stderr, "cmd": " ".join(cmd)})
        # 记录完成/失败日志
        if payload["success"]:
            logger.info(f"STAR完成：samples={len(sample_inputs)} results={results_dir}")
        else:
            logger.warning(f"STAR失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        return payload

    except Exception as e:
        logger.error(f"STAR异常：{e}")
        return {"success": False, "error": f"执行STAR比对失败: {str(e)}"}


@tool
def parse_star_metrics(results_directory: str) -> Dict[str, Any]:
    """解析STAR比对结果文件，提取比对指标
    
    Args:
        results_directory: STAR结果目录路径
    
    Returns:
        Dict: 解析后的比对指标数据
        {
            "success": bool,
            "total_samples": int,
            "results_directory": str,
            "per_sample_metrics": List[Dict],    # 每个样本的详细指标
            "overall_statistics": Dict,          # 总体统计信息
            "quality_assessment": Dict           # 质量评估
        }
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {
                "success": False,
                "error": f"STAR结果目录不存在: {results_directory}"
            }
        
        sample_metrics = []
        overall_stats = {
            "total_input_reads": 0,
            "total_mapped_reads": 0,
            "total_uniquely_mapped": 0,
            "total_multi_mapped": 0,
            "mapping_rates": [],
            "unique_mapping_rates": [],
            "multi_mapping_rates": [],
            "mismatch_rates": []
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
        
        # 遍历样本目录，查找Log.final.out文件
        for sample_dir in star_dir.iterdir():
            if not sample_dir.is_dir():
                continue
                
            # 与 star.nf 中 outFileNamePrefix 保持一致：{sample_id}/{sample_id}.*
            log_final_file = sample_dir / f"{sample_dir.name}.Log.final.out"
            if not log_final_file.exists():
                sample_metrics.append({
                    "sample_id": sample_dir.name,
                    "error": "Log.final.out文件不存在"
                })
                continue
            
            try:
                # 解析Log.final.out文件
                with open(log_final_file, 'r', encoding='utf-8') as f:
                    log_content = f.read()
                
                # 提取关键指标（使用正则表达式）
                input_reads = _extract_metric(log_content, r"Number of input reads.*?(\d+)")
                uniquely_mapped = _extract_metric(log_content, r"Uniquely mapped reads number.*?(\d+)")
                multi_mapped = _extract_metric(log_content, r"Number of reads mapped to multiple loci.*?(\d+)")
                unmapped = _extract_metric(log_content, r"Number of reads unmapped.*?(\d+)")
                
                # 计算比率
                if input_reads > 0:
                    mapping_rate = (uniquely_mapped + multi_mapped) / input_reads
                    unique_mapping_rate = uniquely_mapped / input_reads
                    multi_mapping_rate = multi_mapped / input_reads
                else:
                    mapping_rate = unique_mapping_rate = multi_mapping_rate = 0
                
                # 提取mismatch率
                mismatch_rate = _extract_metric(log_content, r"Mismatch rate per base.*?([\d.]+)%") / 100
                
                sample_metric = {
                    "sample_id": sample_dir.name,
                    "input_reads": input_reads,
                    "uniquely_mapped": uniquely_mapped,
                    "multi_mapped": multi_mapped,
                    "unmapped": unmapped,
                    "mapping_rate": round(mapping_rate, 4),
                    "unique_mapping_rate": round(unique_mapping_rate, 4),
                    "multi_mapping_rate": round(multi_mapping_rate, 4),
                    "mismatch_rate": round(mismatch_rate, 4),
                    "log_file": str(log_final_file)
                }
                
                sample_metrics.append(sample_metric)
                
                # 累加到总体统计
                overall_stats["total_input_reads"] += input_reads
                overall_stats["total_mapped_reads"] += (uniquely_mapped + multi_mapped)
                overall_stats["total_uniquely_mapped"] += uniquely_mapped
                overall_stats["total_multi_mapped"] += multi_mapped
                overall_stats["mapping_rates"].append(mapping_rate)
                overall_stats["unique_mapping_rates"].append(unique_mapping_rate)
                overall_stats["multi_mapping_rates"].append(multi_mapping_rate)
                overall_stats["mismatch_rates"].append(mismatch_rate)
                
            except Exception as e:
                sample_metrics.append({
                    "sample_id": sample_dir.name,
                    "log_file": str(log_final_file),
                    "error": f"解析失败: {str(e)}"
                })
        
        # 计算总体指标
        total_samples = len([m for m in sample_metrics if "error" not in m])
        if total_samples > 0:
            overall_mapping_rate = overall_stats["total_mapped_reads"] / overall_stats["total_input_reads"]
            overall_unique_rate = overall_stats["total_uniquely_mapped"] / overall_stats["total_input_reads"]
            overall_multi_rate = overall_stats["total_multi_mapped"] / overall_stats["total_input_reads"]
            avg_mismatch_rate = sum(overall_stats["mismatch_rates"]) / len(overall_stats["mismatch_rates"])
        else:
            overall_mapping_rate = overall_unique_rate = overall_multi_rate = avg_mismatch_rate = 0
        
        # 质量评估
        quality_assessment = {
            "overall_quality": "good" if overall_mapping_rate > 0.85 else "moderate" if overall_mapping_rate > 0.7 else "poor",
            "unique_mapping_status": "good" if overall_unique_rate > 0.8 else "moderate" if overall_unique_rate > 0.6 else "poor",
            "multi_mapping_status": "good" if overall_multi_rate < 0.2 else "moderate" if overall_multi_rate < 0.3 else "high"
        }
        
        result = {
            "success": True,
            "total_samples": total_samples,
            "results_directory": results_directory,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_input_reads": overall_stats["total_input_reads"],
                "total_mapped_reads": overall_stats["total_mapped_reads"],
                "total_uniquely_mapped": overall_stats["total_uniquely_mapped"],
                "total_multi_mapped": overall_stats["total_multi_mapped"],
                "overall_mapping_rate": round(overall_mapping_rate, 4),
                "overall_unique_mapping_rate": round(overall_unique_rate, 4),
                "overall_multi_mapping_rate": round(overall_multi_rate, 4),
                "average_mismatch_rate": round(avg_mismatch_rate, 4)
            },
            "quality_assessment": quality_assessment
        }
        try:
            logger.info(
                f"STAR结果: samples={result['total_samples']} map={result['overall_statistics']['overall_mapping_rate']} "
                f"unique={result['overall_statistics']['overall_unique_mapping_rate']} multi={result['overall_statistics']['overall_multi_mapping_rate']}"
            )
            if sample_metrics:
                logger.debug(f"STAR样本预览：{sample_metrics[0]}")
        except Exception:
            pass
        return result
        
    except Exception as e:
        try:
            logger.error(f"解析STAR结果失败：{e}")
        except Exception:
            pass
        return {
            "success": False,
            "error": f"解析STAR结果失败: {str(e)}"
        }

def _extract_metric(text: str, pattern: str) -> float:
    """从文本中提取数值指标的辅助函数"""
    import re
    match = re.search(pattern, text)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return 0.0
    return 0.0


@tool
def run_nextflow_hisat2(
    hisat2_params: Dict[str, Any],
    fastp_results: Dict[str, Any],
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None
) -> Dict[str, Any]:
    """执行 HISAT2 比对（精简版，与 run_nextflow_star 等价）

    约束（与路径契约一致）:
    - 仅在 fastp_results.success 为真且包含 per_sample_outputs 时放行
    - 统一复用 FastP 的 results_dir 作为运行根目录
    - HISAT2 索引优先使用 genome_info.hisat2_index_dir；否则由 genome_info.fasta_path 推导
    - sample_inputs 仅来源于 fastp_results.per_sample_outputs（不再扫描目录）
    - per_sample_outputs 路径与 hisat2.nf 产出一致（样本子目录 + 默认文件名）
    """
    try:
        tools_config = get_tools_config()

        # 1) 校验 FastP 结果与运行根目录
        if not (fastp_results and fastp_results.get("success")):
            return {"success": False, "error": "FastP结果无效，无法执行HISAT2比对"}

        fastp_results_dir = fastp_results.get("results_dir")
        if not fastp_results_dir:
            return {"success": False, "error": "FastP结果缺少results_dir"}

        per_sample = fastp_results.get("per_sample_outputs") or []
        if not per_sample:
            return {"success": False, "error": "FastP结果缺少per_sample_outputs"}

        # 2) 运行根目录与工作目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = Path(fastp_results_dir)
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"hisat2_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)

        # 3) 解析 HISAT2 索引目录
        def _resolve(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (tools_config.settings.project_root / pp)

        hisat2_index_prefix = ""
        if isinstance(genome_info, dict):
            hisat2_index_dir = genome_info.get("hisat2_index_dir") or genome_info.get("index_dir") or ""
            if hisat2_index_dir:
                hisat2_index_prefix = str(_resolve(hisat2_index_dir) / "genome")
            else:
                fasta_path = genome_info.get("fasta_path") or genome_info.get("fasta")
                if fasta_path:
                    hisat2_index_prefix = str(tools_config.get_hisat2_index_dir(_resolve(fasta_path)) / "genome")

        if not hisat2_index_prefix:
            return {"success": False, "error": "缺少HISAT2索引目录（genome_info.hisat2_index_dir 或 fasta_path 必须提供）"}

        # 检查索引文件是否存在
        index_files = list(Path(hisat2_index_prefix).parent.glob(f"{Path(hisat2_index_prefix).name}.*.ht2"))
        if not index_files:
            return {"success": False, "error": f"HISAT2索引文件不存在: {hisat2_index_prefix}.*.ht2"}

        # 4) 构造 sample_inputs（仅使用 FastP 返回的结构）
        sample_inputs: List[Dict[str, Any]] = []
        for i, fp in enumerate(per_sample):
            sid = fp.get("sample_id", f"sample_{i+1}")
            r1 = fp.get("trimmed_single") or fp.get("trimmed_r1")
            r2 = fp.get("trimmed_r2")
            if not r1:
                continue
            sample_inputs.append({
                "sample_id": sid,
                "is_paired": bool(r2),
                "read1": r1,
                **({"read2": r2} if r2 else {})
            })
        if not sample_inputs:
            return {"success": False, "error": "未从FastP结果构造到任何样本输入"}

        # 5) 组装 Nextflow 参数
        cleaned_params: Dict[str, Any] = {}
        for k, v in (hisat2_params or {}).items():
            if v is None or k in {"hisat2_cpus", "threads", "p"}:
                continue
            cleaned_params[k.lstrip('-')] = v

        nf_params = {
            "sample_inputs": json.dumps(sample_inputs, ensure_ascii=False),
            "hisat2_index": hisat2_index_prefix,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **cleaned_params,
        }

        # 保存参数文件到hisat2子目录
        hisat2_dir = results_dir / "hisat2"
        hisat2_dir.mkdir(parents=True, exist_ok=True)
        params_file = hisat2_dir / "hisat2_params.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "hisat2.nf",
            Path("/src/nextflow/hisat2.nf"),
        ]
        nextflow_script = next((p for p in nf_candidates if p.exists()), None)
        if nextflow_script is None:
            return {"success": False, "error": "未找到 hisat2.nf 脚本，请检查 /src/nextflow/hisat2.nf 路径", "searched": [str(p) for p in nf_candidates]}

        logger.info(f"执行HISAT2比对 - 参数文件: {params_file}")
        logger.info(f"HISAT2索引: {nf_params['hisat2_index']}")
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, cwd=tools_config.settings.project_root)

        # 7) 组装每样本输出路径（与 hisat2.nf publishDir 对齐）
        hisat2_out = results_dir / "hisat2"
        per_sample_outputs: List[Dict[str, Any]] = []
        for item in sample_inputs:
            sid = item["sample_id"]
            sdir = hisat2_out / sid
            entry = {
                "sample_id": sid,
                "aligned_bam": str(sdir / f"{sid}.hisat2.bam"),
                "align_summary": str(sdir / f"{sid}.align_summary.txt"),
                "bam_index": str(sdir / f"{sid}.hisat2.bam.bai"),
            }
            per_sample_outputs.append(entry)

        payload = {
            "success": result.returncode == 0,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": len(sample_inputs),
            "per_sample_outputs": per_sample_outputs,
        }
        if get_tools_config().settings.debug_mode:
            payload.update({"stdout": result.stdout, "stderr": result.stderr, "cmd": " ".join(cmd)})
        return payload

    except Exception as e:
        return {"success": False, "error": f"执行HISAT2比对失败: {str(e)}"}


@tool
def build_hisat2_index(
    genome_id: str,
    p: Optional[int] = None,
    force_rebuild: bool = False,
    results_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """构建 HISAT2 索引（等价 build_star_index）"""
    try:
        tools_config = get_tools_config()

        # 1) 从基因组配置获取基因组信息
        genome_configs = tools_config.get_genome_configs()
        if genome_id not in genome_configs:
            return {"success": False, "error": f"基因组ID '{genome_id}' 未找到，可用ID: {list(genome_configs.keys())}"}
        
        genome_config = genome_configs[genome_id]
        fasta_path = genome_config.get("fasta_path")
        gtf_path = genome_config.get("gtf_path", "")  # GTF可选
        
        if not fasta_path:
            return {"success": False, "error": f"基因组 '{genome_id}' 配置缺少fasta_path"}

        # 2) 解析路径
        def _resolve(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (tools_config.settings.project_root / pp)

        fasta_file = _resolve(fasta_path)
        gtf_file = _resolve(gtf_path)
        
        # 检查必需文件是否存在
        if not fasta_file.exists():
            return {"success": False, "error": f"FASTA文件不存在: {fasta_file}，请先下载基因组文件"}
        if gtf_file.exists():
            return {"success": False, "error": f"GTF文件不存在: {gtf_file}，请先下载基因组文件"}

        # 3) 确定索引目录
        index_dir = tools_config.get_hisat2_index_dir(fasta_file)
        index_dir.mkdir(parents=True, exist_ok=True)

        # 4) 检查是否需要重建
        index_files = list(index_dir.glob("genome.*.ht2"))
        if index_files and not force_rebuild:
            return {
                "success": True,
                "hisat2_index_dir": str(index_dir),
                "status": "已存在",
                "index_files": [str(f) for f in index_files],
                "message": f"HISAT2索引已存在，跳过构建: {index_dir}"
            }

        # 5) 执行索引构建
        work_dir_name = f"hisat2_index_{genome_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        work_dir = tools_config.settings.data_dir / "work" / work_dir_name
        work_dir.mkdir(parents=True, exist_ok=True)

        nf_params = {
            "genome_fasta": str(fasta_file),
            "genome_gtf": str(_resolve(gtf_path)) if gtf_path else "",
            "hisat2_index_dir": str(index_dir),
            "index_basename": "genome",
            "p": p or 4,
        }

        params_file = work_dir / "build_hisat2_index_params.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "build_hisat2_index.nf",
            Path("/src/nextflow/build_hisat2_index.nf"),
        ]
        nextflow_script = next((p for p in nf_candidates if p.exists()), None)
        if nextflow_script is None:
            return {"success": False, "error": "未找到 build_hisat2_index.nf 脚本"}

        logger.info(f"构建HISAT2索引 - 基因组: {genome_id}")
        logger.info(f"FASTA: {fasta_file}")
        logger.info(f"索引目录: {index_dir}")
        
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400, cwd=tools_config.settings.project_root)

        # 7) 检查结果
        final_index_files = list(index_dir.glob("genome.*.ht2"))
        success = result.returncode == 0 and len(final_index_files) > 0

        payload = {
            "success": success,
            "hisat2_index_dir": str(index_dir),
            "genome_id": genome_id,
            "index_files": [str(f) for f in final_index_files],
            "work_dir": str(work_dir),
            "params_file": str(params_file),
        }

        if success:
            payload["status"] = "构建成功"
            payload["message"] = f"HISAT2索引构建完成: {index_dir}"
        else:
            payload["status"] = "构建失败"
            payload["error"] = f"HISAT2索引构建失败，returncode: {result.returncode}"

        if get_tools_config().settings.debug_mode:
            payload.update({"stdout": result.stdout, "stderr": result.stderr, "cmd": " ".join(cmd)})

        return payload

    except Exception as e:
        return {"success": False, "error": f"构建HISAT2索引失败: {str(e)}"}


@tool
def parse_hisat2_metrics(results_directory: str) -> Dict[str, Any]:
    """解析HISAT2比对结果文件，提取比对指标
    
    Args:
        results_directory: HISAT2结果目录路径
        
    Returns:
        Dict: 包含样本指标和整体统计的结果
    """
    try:
        results_path = Path(results_directory)
        
        if not results_path.exists():
            return {"success": False, "error": f"结果目录不存在: {results_directory}"}
        
        # 查找HISAT2结果子目录
        hisat2_dir = results_path / "hisat2"
        if not hisat2_dir.exists():
            return {"success": False, "error": f"HISAT2子目录不存在: {hisat2_dir}"}
        
        sample_metrics = []
        overall_stats = {
            "total_input_reads": 0,
            "total_mapped_reads": 0,
            "total_uniquely_mapped": 0,
            "total_multi_mapped": 0,
            "mapping_rates": [],
            "unique_mapping_rates": [],
            "multi_mapping_rates": []
        }
        
        # 扫描样本目录
        sample_dirs = [d for d in hisat2_dir.iterdir() if d.is_dir()]
        total_samples = len(sample_dirs)
        
        if total_samples == 0:
            return {"success": False, "error": "未找到任何样本目录"}
        
        for sample_dir in sample_dirs:
            sample_id = sample_dir.name
            summary_file = sample_dir / f"{sample_id}.align_summary.txt"
            
            if not summary_file.exists():
                continue
            
            # 解析HISAT2比对摘要
            with open(summary_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # HISAT2输出格式解析
            metrics = _parse_hisat2_summary(content, sample_id)
            if metrics:
                sample_metrics.append(metrics)
                
                # 累积整体统计
                overall_stats["total_input_reads"] += metrics["input_reads"]
                overall_stats["total_mapped_reads"] += metrics["mapped_reads"]
                overall_stats["total_uniquely_mapped"] += metrics.get("uniquely_mapped", 0)
                overall_stats["total_multi_mapped"] += metrics.get("multi_mapped", 0)
                overall_stats["mapping_rates"].append(metrics["mapping_rate"])
                overall_stats["unique_mapping_rates"].append(metrics.get("unique_mapping_rate", 0))
                overall_stats["multi_mapping_rates"].append(metrics.get("multi_mapping_rate", 0))
        
        # 计算整体比率
        if overall_stats["total_input_reads"] > 0:
            overall_mapping_rate = overall_stats["total_mapped_reads"] / overall_stats["total_input_reads"]
            overall_unique_rate = overall_stats["total_uniquely_mapped"] / overall_stats["total_input_reads"]
            overall_multi_rate = overall_stats["total_multi_mapped"] / overall_stats["total_input_reads"]
        else:
            overall_mapping_rate = overall_unique_rate = overall_multi_rate = 0
        
        # 质量评估
        quality_assessment = {
            "overall_quality": "good" if overall_mapping_rate > 0.85 else "moderate" if overall_mapping_rate > 0.7 else "poor",
            "unique_mapping_status": "good" if overall_unique_rate > 0.8 else "moderate" if overall_unique_rate > 0.6 else "poor",
            "multi_mapping_status": "good" if overall_multi_rate < 0.2 else "moderate" if overall_multi_rate < 0.3 else "high"
        }
        
        return {
            "success": True,
            "total_samples": total_samples,
            "results_directory": results_directory,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_input_reads": overall_stats["total_input_reads"],
                "total_mapped_reads": overall_stats["total_mapped_reads"],
                "total_uniquely_mapped": overall_stats["total_uniquely_mapped"],
                "total_multi_mapped": overall_stats["total_multi_mapped"],
                "overall_mapping_rate": round(overall_mapping_rate, 4),
                "overall_unique_mapping_rate": round(overall_unique_rate, 4),
                "overall_multi_mapping_rate": round(overall_multi_rate, 4),
            },
            "quality_assessment": quality_assessment
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"解析HISAT2结果失败: {str(e)}"
        }


def _parse_hisat2_summary(content: str, sample_id: str) -> Dict[str, Any]:
    """解析HISAT2比对摘要内容（健壮版）

    兼容数字包含千分位逗号、单端/双端输出差异，支持从 overall alignment rate 回退计算。
    """
    try:
        lines = [ln.strip() for ln in (content or "").splitlines() if ln.strip()]
        metrics: Dict[str, Any] = {"sample_id": sample_id}

        # 工具函数：数字提取（支持1,234,567）
        def _parse_int_token(tok: str) -> int:
            try:
                return int(str(tok).replace(",", "").strip())
            except Exception:
                return 0

        # 提取总reads
        # 例："25,000,000 reads; of these:" 或 "25000000 reads; of these:"
        for ln in lines:
            if "reads; of these:" in ln:
                # 优先用正则
                m = re.search(r"([\d,]+)\s+reads; of these:", ln)
                if m:
                    metrics["input_reads"] = _parse_int_token(m.group(1))
                else:
                    # 回退：取首个token
                    parts = ln.split()
                    if parts:
                        metrics["input_reads"] = _parse_int_token(parts[0])
                break

        # 提取映射计数（双端优先匹配concordantly，单端匹配非concordantly）
        for ln in lines:
            if "aligned concordantly exactly 1 time" in ln:
                metrics["uniquely_mapped"] = _extract_reads_count(ln)
            elif "aligned concordantly >1 times" in ln:
                metrics["multi_mapped"] = _extract_reads_count(ln)
            elif ("aligned exactly 1 time" in ln) and ("concordantly" not in ln):
                metrics["uniquely_mapped"] = _extract_reads_count(ln)
            elif ("aligned >1 times" in ln) and ("concordantly" not in ln):
                metrics["multi_mapped"] = _extract_reads_count(ln)

        # 解析 overall alignment rate（百分比），作为映射率的兜底
        overall_rate = None
        for ln in lines[::-1]:
            # 例："90.00% overall alignment rate"
            if "overall alignment rate" in ln:
                m = re.search(r"([\d.]+)%\s+overall alignment rate", ln)
                if m:
                    try:
                        overall_rate = float(m.group(1)) / 100.0
                    except Exception:
                        overall_rate = None
                break

        # 计算派生指标
        if "input_reads" in metrics:
            total_reads = int(metrics.get("input_reads", 0) or 0)
            uniquely_mapped = int(metrics.get("uniquely_mapped", 0) or 0)
            multi_mapped = int(metrics.get("multi_mapped", 0) or 0)
            mapped_reads = uniquely_mapped + multi_mapped

            # 若unique+multi均未捕获，但有overall_rate，回退计算mapped_reads
            if mapped_reads == 0 and overall_rate is not None and total_reads > 0:
                mapped_reads = int(round(total_reads * overall_rate))

            unmapped_reads = max(0, total_reads - mapped_reads)
            metrics["mapped_reads"] = mapped_reads
            metrics["unmapped_reads"] = unmapped_reads

            if total_reads > 0:
                metrics["mapping_rate"] = round(mapped_reads / total_reads, 4)
                metrics["unique_mapping_rate"] = round(uniquely_mapped / total_reads, 4) if uniquely_mapped else 0.0
                metrics["multi_mapping_rate"] = round(multi_mapped / total_reads, 4) if multi_mapped else 0.0
                metrics["unmapped_rate"] = round(unmapped_reads / total_reads, 4)
            else:
                metrics.update({
                    "mapping_rate": 0.0,
                    "unique_mapping_rate": 0.0,
                    "multi_mapping_rate": 0.0,
                    "unmapped_rate": 0.0,
                })

        return metrics if "input_reads" in metrics else None

    except Exception as e:
        logger.error(f"解析HISAT2摘要失败 (样本 {sample_id}): {e}")
        return None


def _extract_reads_count(line: str) -> int:
    """从HISAT2输出行中提取reads数量（支持千分位逗号）"""
    match = re.search(r'([\d,]+)\s+\([\d.]+%\)', line.strip())
    if match:
        try:
            return int(match.group(1).replace(',', ''))
        except Exception:
            return 0
    return 0


# ==================== FeatureCounts 专用工具函数 ====================

@tool
def run_nextflow_featurecounts(
    featurecounts_params: Dict[str, Any],
    star_results: Dict[str, Any],
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None,
    base_results_dir: Optional[str] = None,
    hisat2_results: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """执行Nextflow FeatureCounts定量流程

    Args:
        featurecounts_params: FeatureCounts参数字典（可包含 -T/-s/-p/-M 等风格键）
        star_results: STAR节点结果（可为空），需包含 per_sample_outputs 中的 BAM 路径
        genome_info: 基因组信息（需提供 GTF 路径，如 gtf_path）
        results_timestamp: 可选的时间戳，优先用于结果目录
        base_results_dir: 可选的基底结果目录（来自Detect节点）
        hisat2_results: HISAT2节点结果（可为空），与 star_results 二选一

    Returns:
        执行结果字典，包含状态、输出路径、样本输出等
    """
    try:
        tools_config = get_tools_config()

        # 仅支持容器路径的归一化（不做本地映射）
        def _to_container_path(path_str: str) -> str:
            """将输入规范化为容器内路径：
            - 允许 '/data/...', '/config/...', '/work/...'
            - 'genomes/...' 将映射为 '/data/genomes/...'
            - 其他相对或本地路径一律原样返回（由存在性校验报错）
            """
            s = str(path_str or '').strip()
            if not s:
                return s
            if s.startswith('/data/') or s.startswith('/config/') or s.startswith('/work/'):
                return s
            if s.startswith('genomes/'):
                return '/data/' + s
            return s

        # 校验依赖输入
        # 选择可用的比对结果（STAR/HISAT2）
        align_results = None
        if star_results and star_results.get("success"):
            align_results = star_results
        elif hisat2_results and hisat2_results.get("success"):
            align_results = hisat2_results
        else:
            return {"success": False, "error": "比对结果无效（STAR/HISAT2），无法执行FeatureCounts"}

        per_sample = align_results.get("per_sample_outputs") or []
        if not per_sample:
            return {"success": False, "error": "比对结果缺少 per_sample_outputs 信息"}

        # 环境检查：允许本地与容器环境，路径规范化在下方处理

        # 解析并归一化 GTF 注释文件路径（容器内）
        gtf_file_raw = (
            genome_info.get("gtf_path")
            or genome_info.get("gtf")
            or genome_info.get("annotation_gtf")
            or ""
        )
        if not gtf_file_raw:
            return {"success": False, "error": "genome_info 未提供 GTF 注释文件路径 (gtf_path)"}
        gtf_file = _to_container_path(gtf_file_raw)
        if not Path(gtf_file).exists():
            return {
                "success": False,
                "error": f"GTF文件不存在: {gtf_file}",
            }

        # 运行根目录（results_dir）：复用比对步骤的 results_dir，保持同一运行根目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        run_root = Path(align_results.get("results_dir") or base_results_dir or tools_config.results_dir / f"{timestamp}")
        # 容器路径规范（如需）
        run_root = Path(_to_container_path(str(run_root)))
        results_dir = run_root
        # 统一 Nextflow 工作目录到 /data/work，使用运行ID区分
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"featurecounts_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        (results_dir / "featurecounts").mkdir(parents=True, exist_ok=True)
        logger.info(f"FeatureCounts启动：bam={len(per_sample)} results={results_dir} work={work_dir}")

        # 构建 Nextflow 参数（与 featurecounts.nf 对齐）
        # 将 STAR 输出的 BAM 列表转换为JSON字符串（Nextflow 端会 echo 后再解析）
        bam_entries = []
        for item in per_sample:
            sid = item.get("sample_id") or "sample"
            bam = item.get("aligned_bam") or item.get("bam")
            if not bam:
                continue
            bam_norm = _to_container_path(bam)
            if not Path(bam_norm).exists():
                return {
                    "success": False,
                    "error": f"BAM文件不存在: {bam_norm}",
                }
            bam_entries.append({"sample_id": sid, "bam_file": bam_norm})
        if not bam_entries:
            return {"success": False, "error": "未从STAR结果中收集到任何BAM路径"}

        # 参数映射：Python风格/短旗标 → Nextflow params 名称
        mapped: Dict[str, Any] = {}
        p = featurecounts_params or {}

        def pick_bool(key: str) -> Optional[bool]:
            v = p.get(key)
            if isinstance(v, bool):
                return v
            return None

        def pick_int(key: str) -> Optional[int]:
            v = p.get(key)
            try:
                return int(v) if v is not None else None
            except Exception:
                return None

        # 线程/链特异性/特征/属性/质量阈
        mapped["threads"] = pick_int("-T") or p.get("threads") or 4
        mapped["strand_specificity"] = pick_int("-s") or p.get("strand_specificity") or 0
        mapped["feature_type"] = p.get("-t") or p.get("feature_type") or "exon"
        mapped["attribute_type"] = p.get("-g") or p.get("attribute_type") or "gene_id"
        mapped["min_mapping_quality"] = pick_int("-Q") or p.get("min_mapping_quality") or 10

        # 布尔开关 - 修改count_reads_pairs默认值为false
        mapped["count_reads_pairs"] = pick_bool("-p") if pick_bool("-p") is not None else (p.get("count_reads_pairs") if isinstance(p.get("count_reads_pairs"), bool) else False)
        mapped["count_multi_mapping_reads"] = pick_bool("-M") if pick_bool("-M") is not None else bool(p.get("count_multi_mapping_reads", False))
        mapped["ignore_duplicates"] = bool(p.get("--ignoreDup", False) or p.get("ignore_duplicates", False))
        mapped["require_both_ends_mapped"] = bool(p.get("-B", False) or p.get("require_both_ends_mapped", False))
        mapped["exclude_chimeric"] = bool(p.get("-C", False) or p.get("exclude_chimeric", False))

        # 组装 Nextflow 参数文件
        nf_params = {
            "input_bam_list": json.dumps(bam_entries, ensure_ascii=False),
            "gtf_file": gtf_file,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **mapped,
        }

        # 保存参数文件到featurecounts子目录
        fc_dir = results_dir / "featurecounts"
        fc_dir.mkdir(parents=True, exist_ok=True)
        params_file = fc_dir / "featurecounts_params.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 定位 Nextflow 脚本
        nf_candidates = [
            tools_config.settings.project_root / "src" / "nextflow" / "featurecounts.nf",
            Path("/src/nextflow/featurecounts.nf"),
        ]
        nextflow_script = None
        for cand in nf_candidates:
            if cand.exists():
                nextflow_script = cand
                break
        if nextflow_script is None:
            return {
                "success": False,
                "error": "未找到 featurecounts.nf 脚本，请检查 /src/nextflow/featurecounts.nf 路径",
                "searched": [str(p) for p in nf_candidates],
            }

        # 执行 Nextflow
        cmd = [
            "nextflow",
            "run",
            str(nextflow_script),
            "-params-file",
            str(params_file),
            "-work-dir",
            str(work_dir),
        ]

        logger.info("执行Nextflow FeatureCounts流水线")
        logger.info(f"参数文件: {params_file}")
        logger.info(f"工作目录: {work_dir}")
        logger.info(f"结果目录: {results_dir}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,
            cwd=tools_config.settings.project_root,
        )

        # 构建输出结构 - 适配新的批量输出格式
        sample_count = len(bam_entries)
        per_sample_outputs: List[Dict[str, Any]] = []
        fc_root = results_dir / "featurecounts"
        
        # 新的featurecounts.nf脚本生成批量文件，不再有每个样本的单独目录
        # 主要输出文件：
        # - all_samples.featureCounts (完整计数矩阵)
        # - all_samples.featureCounts.summary (统计摘要)
        # - merged_counts_matrix.txt (兼容格式的矩阵)
        
        # 为兼容性生成per_sample_outputs结构，指向批量文件
        for entry in bam_entries:
            sid = entry["sample_id"]
            sample_output = {
                "sample_id": sid,
                "counts_file": str(fc_root / "all_samples.featureCounts"),  # 指向批量文件
                "summary_file": str(fc_root / "all_samples.featureCounts.summary"),  # 指向批量文件
            }
            per_sample_outputs.append(sample_output)

        payload = {
            "success": result.returncode == 0,
            "message": f"FeatureCounts定量完成，处理了{sample_count}个样本" if result.returncode == 0 else "FeatureCounts执行失败",
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": sample_count,
            "per_sample_outputs": per_sample_outputs,
            "matrix_path": str(fc_root / "merged_counts_matrix.txt"),
            "nextflow_params": nf_params,
        }
        try:
            if get_tools_config().settings.debug_mode:
                payload.update({
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                    "cmd": " ".join(cmd),
                })
            else:
                # 非调试模式去掉 nextflow_params 以减小负载
                payload.pop("nextflow_params", None)
        except Exception:
            pass
        if payload["success"]:
            logger.info(f"FeatureCounts完成：samples={sample_count} results={results_dir}")
        else:
            logger.warning(f"FeatureCounts失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        return payload

    except subprocess.TimeoutExpired:
        logger.warning("FeatureCounts执行超时：30分钟")
        return {
            "success": False,
            "error": "Nextflow执行超时（30分钟）",
        }
    except Exception as e:
        logger.error(f"FeatureCounts异常：{e}")
        return {
            "success": False,
            "error": f"执行FeatureCounts流水线失败: {str(e)}",
        }


@tool
def parse_featurecounts_metrics(results_directory: str) -> Dict[str, Any]:
    """解析FeatureCounts定量结果，输出样本级与总体指标
    
    Args:
        results_directory: FeatureCounts结果目录（包含 featurecounts 子目录）
    
    Returns:
        解析后的指标（assignment rates、未分配原因等）
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {"success": False, "error": f"结果目录不存在: {results_directory}"}
        
        fc_dir = results_path / "featurecounts"
        if not fc_dir.exists():
            return {"success": False, "error": f"缺少特征计数目录: {fc_dir}"}
        
        # 查找批量输出的汇总文件
        summary_file = fc_dir / "all_samples.featureCounts.summary"
        counts_file = fc_dir / "all_samples.featureCounts"
        
        if not summary_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts汇总文件: {summary_file}"}
        
        if not counts_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts计数文件: {counts_file}"}
        
        # 解析批量汇总文件
        sample_metrics: List[Dict[str, Any]] = []
        totals = {
            "assigned": 0,
            "nofeatures": 0,
            "multimapping": 0,
            "ambiguous": 0,
            "mappingquality": 0,
            "other": 0,
            "total": 0,
        }
        
        try:
            # 规范化样本ID的内部工具：
            # - 兼容列名为BAM文件路径/文件名/带STAR后缀的多种情况
            # - 目标：与样本ID（如 SRRxxxx、样本目录名）对齐
            def _normalize_sample_id(name: str) -> str:
                s = str(name or "").strip()
                if not s:
                    return s
                # 去除可能的路径前缀（同时支持 / 与 \\ 分隔符）
                if "/" in s:
                    s = s.split("/")[-1]
                if "\\" in s:
                    s = s.split("\\")[-1]
                # 去除常见扩展名
                for ext in [".bam", ".cram", ".sam", ".txt"]:
                    if s.endswith(ext):
                        s = s[: -len(ext)]
                        break
                # 去除常见后缀（STAR/HISAT2 的命名后缀，点/下划线变体）
                star_suffixes = [
                    ".Aligned.sortedByCoord.out",
                    ".Aligned.out",
                    ".Aligned",
                    "_Aligned.sortedByCoord.out",
                    "_Aligned.out",
                    "_Aligned",
                    ".hisat2",  # 来自 HISAT2 的常见后缀（在移除 .bam 后可能残留）
                    "_hisat2",
                ]
                for suf in star_suffixes:
                    if s.endswith(suf):
                        s = s[: -len(suf)]
                        break
                return s

            with open(summary_file, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
                
                if len(lines) < 2:
                    return {"success": False, "error": "汇总文件格式错误"}
                
                # 解析标题行获取样本名称（用于回退）
                header = lines[0].strip().split('\t')
                if len(header) < 2:
                    return {"success": False, "error": "汇总文件标题行格式错误"}

                # 优先从参数文件读取样本ID顺序（与执行输入一致）
                sample_ids: List[str] = []
                try:
                    params_path = results_path / "featurecounts" / "featurecounts_params.json"
                    if params_path.exists():
                        with open(params_path, "r", encoding="utf-8") as pf:
                            pf_json = json.load(pf)
                        input_bam_list = pf_json.get("input_bam_list")
                        if isinstance(input_bam_list, str):
                            input_bam_list = json.loads(input_bam_list)
                        if isinstance(input_bam_list, list):
                            for ent in input_bam_list:
                                sid = ent.get("sample_id") or _normalize_sample_id(ent.get("bam_file", ""))
                                sample_ids.append(_normalize_sample_id(sid))
                except Exception:
                    # 若读取或解析失败，忽略并回退到header
                    sample_ids = []

                # 校验样本数是否与汇总列数一致；否则回退到header列名
                if not sample_ids or len(sample_ids) != (len(header) - 1):
                    sample_names = header[1:]  # 第一列是Status，后面是样本名（通常为输入BAM的文件名）
                    sample_ids = [_normalize_sample_id(nm) for nm in sample_names]

                # 初始化每个样本的指标
                for sid in sample_ids:
                    
                    sample_metrics.append({
                        "sample_id": sid,
                        "assigned": 0,
                        "unassigned_unmapped": 0,
                        "unassigned_mappingquality": 0,
                        "unassigned_nofeatures": 0,
                        "unassigned_ambiguity": 0,
                        "total_reads": 0
                    })
                
                # 解析每一行统计数据
                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                        
                    status = parts[0]
                    values = [int(v) for v in parts[1:]]
                    
                    # 更新每个样本的指标
                    for i, value in enumerate(values):
                        if i < len(sample_metrics):
                            if status == "Assigned":
                                sample_metrics[i]["assigned"] = value
                                totals["assigned"] += value
                            elif status == "Unassigned_Unmapped":
                                sample_metrics[i]["unassigned_unmapped"] = value
                            elif status == "Unassigned_MappingQuality":
                                sample_metrics[i]["unassigned_mappingquality"] = value
                                totals["mappingquality"] += value
                            elif status == "Unassigned_NoFeatures":
                                sample_metrics[i]["unassigned_nofeatures"] = value
                                totals["nofeatures"] += value
                            elif status == "Unassigned_Ambiguity":
                                sample_metrics[i]["unassigned_ambiguity"] = value
                                totals["ambiguous"] += value
                
                # 计算每个样本的总读数和分配率
                for sample_metric in sample_metrics:
                    sample_metric["total_reads"] = (
                        sample_metric["assigned"] +
                        sample_metric["unassigned_unmapped"] +
                        sample_metric["unassigned_mappingquality"] +
                        sample_metric["unassigned_nofeatures"] +
                        sample_metric["unassigned_ambiguity"]
                    )
                    
                    if sample_metric["total_reads"] > 0:
                        sample_metric["assignment_rate"] = round(
                            sample_metric["assigned"] / sample_metric["total_reads"], 4
                        )
                    else:
                        sample_metric["assignment_rate"] = 0.0
                    
                    totals["total"] += sample_metric["total_reads"]
        
        except Exception as e:
            return {"success": False, "error": f"解析汇总文件失败: {str(e)}"}
        
        # 计算总体统计
        total_assignment_rate = totals["assigned"] / totals["total"] if totals["total"] > 0 else 0.0
        
        # 读取计数矩阵获取基因数量
        gene_count = 0
        try:
            with open(counts_file, "r", encoding="utf-8") as f:
                # 跳过注释行
                for line in f:
                    if not line.startswith('#'):
                        gene_count += 1
                gene_count -= 1  # 减去标题行
        except Exception:
            gene_count = 0
        
        # 重要文件路径（若存在合并矩阵，则优先提供）
        matrix_file = fc_dir / "merged_counts_matrix.txt"

        return {
            "success": True,
            "results_directory": results_directory,
            "summary_file": str(summary_file),
            "counts_file": str(counts_file),
            "matrix_path": str(matrix_file) if matrix_file.exists() else str(counts_file),
            "sample_count": len(sample_metrics),
            "gene_count": gene_count,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_reads": totals["total"],
                "total_assigned": totals["assigned"],
                "total_unassigned_nofeatures": totals["nofeatures"],
                "total_unassigned_ambiguity": totals["ambiguous"],
                "total_unassigned_mappingquality": totals["mappingquality"],
                "overall_assignment_rate": round(total_assignment_rate, 4)
            }
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": f"解析FeatureCounts结果失败: {str(e)}"
        }


# ==================== Analysis 报告生成工具函数 ====================

@tool
def write_analysis_markdown(
    analysis_report: Dict[str, Any],
    output_dir: Optional[str] = None,
    filename: Optional[str] = None,
    append_llm_section: bool = True
) -> Dict[str, Any]:
    """将结构化的 analysis_report 渲染为可读的 Markdown 摘要文件
    
    Args:
        analysis_report: 按设计文档第4节JSON结构约定生成的分析报告
        output_dir: 目标目录，若为空则使用 context.results_dir/reports/timestamp
        filename: 文件名，默认 analysis_summary.md
        append_llm_section: 是否追加LLM洞察到文末
    
    Returns:
        Dict: 包含成功状态、文件路径、字节数等信息
    """
    try:
        from .config import get_tools_config
        tools_config = get_tools_config()
        
        # 参数验证
        if not analysis_report or not isinstance(analysis_report, dict):
            return {
                "success": False,
                "error": "analysis_report 参数无效或为空"
            }
        
        # 确定输出目录和文件名
        filename = filename or "analysis_summary.md"
        
        if output_dir:
            target_dir = Path(output_dir)
        else:
            # 从report中获取目录信息
            context = analysis_report.get("context", {})
            results_dir = context.get("results_dir")
            timestamp = context.get("timestamp")
            
            if results_dir and timestamp:
                target_dir = Path(results_dir) / "reports" / timestamp
            else:
                # fallback到默认报告目录
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                target_dir = tools_config.reports_dir / timestamp
        
        # 创建目录
        target_dir.mkdir(parents=True, exist_ok=True)
        output_path = target_dir / filename
        
        # 渲染Markdown内容
        markdown_content = _render_analysis_markdown(analysis_report, append_llm_section)
        
        # 写入文件
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(markdown_content)
        
        # 计算统计信息
        file_size = output_path.stat().st_size
        line_count = markdown_content.count('\n') + 1
        
        logger.info(f"写入Markdown报告：{output_path} bytes={file_size} lines={line_count}")
        return {
            "success": True,
            "path": str(output_path.absolute()),
            "bytes": file_size,
            "lines": line_count,
            "used_output_dir": str(target_dir),
            "used_filename": filename
        }
        
    except Exception as e:
        logger.error(f"生成Markdown报告失败：{e}")
        return {
            "success": False,
            "error": f"生成Markdown报告失败: {str(e)}"
        }


def _render_analysis_markdown(report: Dict[str, Any], append_llm: bool = True) -> str:
    """渲染分析报告为Markdown格式的辅助函数"""
    
    def _safe_get(obj: Dict, *keys, default="-"):
        """安全获取嵌套字典值"""
        for key in keys:
            if isinstance(obj, dict) and key in obj:
                obj = obj[key]
            else:
                return default
        return obj if obj is not None else default
    
    def _format_percent(value, default="-"):
        """格式化百分比"""
        if value is None or value == default:
            return default
        try:
            return f"{float(value) * 100:.1f}%" if 0 <= float(value) <= 1 else f"{float(value):.1f}%"
        except (ValueError, TypeError):
            return default
    
    def _format_number(value, default="-"):
        """格式化数字"""
        if value is None:
            return default
        try:
            if isinstance(value, float):
                return f"{value:,.1f}" if value >= 1000 else f"{value:.3f}"
            else:
                return f"{int(value):,}"
        except (ValueError, TypeError):
            return default
    
    # 获取报告各部分
    pipeline = report.get("pipeline", {})
    context = report.get("context", {})
    metrics = report.get("metrics", {})
    per_sample = report.get("per_sample", [])
    summary = report.get("summary", {})
    files = report.get("files", {})
    recommendations = report.get("recommendations", [])
    llm_output = report.get("llm", report.get("llm_output", {}))
    
    # 开始构建Markdown
    lines = []
    
    # 1) 标题与基本信息
    steps_str = " → ".join(pipeline.get("steps", []))
    lines.extend([
        f"# RNA-seq 分析报告",
        "",
        f"**分析流程**: {steps_str}",
        f"**物种**: {pipeline.get('species', '未知')}",
        f"**基因组版本**: {pipeline.get('genome_version', '未知')}",
        f"**样本数量**: {context.get('sample_count', 0)}",
        f"**分析时间**: {context.get('timestamp', '未知')}",
        f"**结果目录**: `{context.get('results_dir', '未知')}`",
        ""
    ])
    
    # 2) 总体结论与关键发现
    status = summary.get("status", "UNKNOWN")
    status_emoji = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}.get(status, "❓")
    
    samples_info = summary.get("samples", {})
    lines.extend([
        f"## {status_emoji} 总体结论: {status}",
        "",
        f"**样本统计**:",
        f"- 通过: {samples_info.get('pass', 0)} 个",
        f"- 警告: {samples_info.get('warn', 0)} 个", 
        f"- 失败: {samples_info.get('fail', 0)} 个",
        ""
    ])
    
    key_findings = summary.get("key_findings", [])
    if key_findings:
        lines.extend([
            f"**关键发现**:",
            *[f"- {finding}" for finding in key_findings],
            ""
        ])
    
    # 3) 各步骤总体指标
    lines.append("## 📊 流程步骤指标")
    lines.append("")
    
    # FastP 指标
    fastp_metrics = metrics.get("fastp", {}).get("overall", {})
    if fastp_metrics:
        lines.extend([
            "### FastP 质量控制",
            f"- 平均Q30质量: {_format_percent(fastp_metrics.get('average_q30_rate'))}",
            f"- 读数通过率: {_format_percent(fastp_metrics.get('overall_read_pass_rate'))}",
            f"- 碱基通过率: {_format_percent(fastp_metrics.get('overall_base_pass_rate'))}",
            f"- 总reads处理: {_format_number(fastp_metrics.get('total_reads_before'))} → {_format_number(fastp_metrics.get('total_reads_after'))}",
            ""
        ])
    
    # STAR 指标
    star_metrics = metrics.get("star", {}).get("overall", {})
    if star_metrics:
        lines.extend([
            "### STAR 序列比对",
            f"- 总体比对率: {_format_percent(star_metrics.get('overall_mapping_rate'))}",
            f"- 唯一比对率: {_format_percent(star_metrics.get('overall_unique_mapping_rate'))}",
            f"- 多重比对率: {_format_percent(star_metrics.get('overall_multi_mapping_rate'))}",
            f"- 平均错配率: {_format_percent(star_metrics.get('average_mismatch_rate'))}",
            ""
        ])
    
    # FeatureCounts 指标  
    fc_metrics = metrics.get("featurecounts", {}).get("overall", {})
    if fc_metrics:
        lines.extend([
            "### FeatureCounts 基因定量",
            f"- 整体分配率: {_format_percent(fc_metrics.get('overall_assignment_rate'))}",
            f"- 总分配reads: {_format_number(fc_metrics.get('total_assigned'))}",
            f"- 未分配(无特征): {_format_number(fc_metrics.get('total_unassigned_nofeatures'))}",
            f"- 未分配(歧义): {_format_number(fc_metrics.get('total_unassigned_ambiguity'))}",
            ""
        ])
    
    # 4) 样本级摘要 - 按健康度排序
    if per_sample:
        lines.extend([
            "## 🔬 样本详情",
            "",
            "| 样本ID | 健康度 | Q30 | 比对率 | 分配率 | 备注 |",
            "|--------|--------|-----|--------|--------|------|"
        ])
        
        # 按健康度排序 (FAIL > WARN > PASS)
        health_order = {"FAIL": 0, "WARN": 1, "PASS": 2}
        sorted_samples = sorted(per_sample, key=lambda x: health_order.get(x.get("health", "UNKNOWN"), 3))
        
        for sample in sorted_samples:
            sid = sample.get("sample_id", "")
            health = sample.get("health", "")
            health_emoji = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}.get(health, "❓")
            
            # 提取关键指标
            fastp_data = sample.get("fastp", {})
            star_data = sample.get("star", {})
            fc_data = sample.get("featurecounts", {})
            
            q30 = _format_percent(fastp_data.get("q30_after"))
            mapping = _format_percent(star_data.get("mapping_rate"))
            assignment = _format_percent(fc_data.get("assignment_rate"))
            
            notes = sample.get("notes", [])
            notes_str = ", ".join(notes[:2]) if notes else "-"
            if len(notes) > 2:
                notes_str += "..."
                
            lines.append(f"| {sid} | {health_emoji} {health} | {q30} | {mapping} | {assignment} | {notes_str} |")
        
        lines.append("")
    
    # 5) 重要文件
    matrix_path = files.get("matrix_path") or _safe_get(fc_metrics, "matrix_path")
    report_json = files.get("report_json")
    
    if matrix_path or report_json:
        lines.extend([
            "## 📁 重要文件",
            ""
        ])
        if matrix_path:
            lines.append(f"- **计数矩阵**: `{matrix_path}`")
        if report_json:
            lines.append(f"- **详细报告**: `{report_json}`")
        lines.append("")
    
    # 6) 建议与后续步骤
    if recommendations:
        lines.extend([
            "## 💡 建议与后续步骤",
            ""
        ])
        
        for rec in recommendations:
            rec_type = rec.get("type", "")
            title = rec.get("title", "")
            detail = rec.get("detail", "")
            
            if rec_type == "action":
                lines.append(f"### 🔧 {title}")
            elif rec_type == "next":
                lines.append(f"### 📈 {title}")
            else:
                lines.append(f"### {title}")
            
            lines.extend([f"{detail}", ""])
    
    # 7) LLM 洞察 (可选)
    if append_llm and llm_output:
        lines.extend([
            "## 🤖 智能分析洞察",
            ""
        ])
        
        global_summary = llm_output.get("global_summary")
        if global_summary:
            lines.extend([
                "### 总体评估",
                global_summary,
                ""
            ])
        
        llm_findings = llm_output.get("key_findings", [])
        if llm_findings:
            lines.extend([
                "### 关键发现",
                *[f"- {finding}" for finding in llm_findings],
                ""
            ])
        
        per_sample_flags = llm_output.get("per_sample_flags", [])
        if per_sample_flags:
            lines.extend([
                "### 样本级问题",
                ""
            ])
            for flag in per_sample_flags:
                sid = flag.get("sample_id", "")
                issues = flag.get("issues", [])
                severity = flag.get("severity", "")
                if issues:
                    lines.append(f"**{sid}** ({severity}): " + ", ".join(issues))
            lines.append("")
        
        risks = llm_output.get("risks", [])
        if risks:
            lines.extend([
                "### ⚠️ 潜在风险",
                *[f"- {risk}" for risk in risks],
                ""
            ])
        
        # 如果有额外的LLM生成的Markdown片段
        report_md = llm_output.get("report_md")
        if report_md:
            lines.extend([
                "### 补充分析",
                report_md,
                ""
            ])
    
    # 生成时间戳
    lines.extend([
        "---",
        f"*报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*"
    ])
    
    return "\n".join(lines)
