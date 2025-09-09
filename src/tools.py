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
    
    return {
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
        
        return {
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

        # 结果目录：优先使用 sample_info 提供；否则按时间戳生成到 data/results 下
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = sample_info.get("results_dir") or str(config.settings.data_dir / "results" / f"fastp_{timestamp}")

        # 工作目录：固定放到 data/tmp 下（包含时间戳/结果名）
        run_id = results_dir.split('/')[-1] if '/' in str(results_dir) else timestamp
        temp_dir = Path(base_data_path) / "tmp" / f"nextflow_fastp_{run_id}"
        temp_dir.mkdir(parents=True, exist_ok=True)
        results_dir = Path(results_dir)
        print(f"📁 运行目录: base={base_data_path} work={temp_dir} results={results_dir}")
        
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
        
        # 创建Nextflow参数文件
        params_file = temp_dir / "fastp_params.json"
        with open(params_file, 'w', encoding='utf-8') as f:
            json.dump(nextflow_params, f, indent=2, ensure_ascii=False)
        
        # 构建Nextflow命令（兼容Docker与本地路径）
        # 1) 优先使用项目根目录下的 fastp.nf（本地开发）
        # 2) Docker镜像中 fastp.nf 位于根路径 '/'（见 Dockerfile COPY fastp.nf /）
        nf_candidates = [
            config.settings.project_root / "fastp.nf",
            Path("/fastp.nf")
        ]
        nextflow_script = None
        for cand in nf_candidates:
            if cand.exists():
                nextflow_script = cand
                break
        if nextflow_script is None:
            return {
                "success": False,
                "error": "未找到 fastp.nf 脚本，请检查容器内是否存在 /fastp.nf 或本地项目根目录",
                "searched": [str(p) for p in nf_candidates]
            }

        cmd = [
            "nextflow", "run",
            str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(temp_dir / "work"),
            "--data", str(config.settings.data_dir)
        ]
        
        print(f"🚀 执行Nextflow FastP流水线...")
        print(f"   参数文件: {params_file}")
        print(f"   工作目录: {temp_dir / 'work'}")
        print(f"   结果目录: {results_dir}")
        
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
            sample_count = len(sample_info["sample_groups"])
            
            return {
                "success": True,
                "message": f"FastP质控完成，处理了{sample_count}个样本",
                "execution_time": execution_time,
                "results_dir": str(results_dir),
                "work_dir": str(temp_dir / "work"),
                "params_file": str(params_file),
                "sample_count": sample_count,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "nextflow_params": nextflow_params
            }
        else:
            return {
                "success": False,
                "error": f"Nextflow执行失败 (返回码: {result.returncode})",
                "execution_time": execution_time,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "cmd": " ".join(cmd)
            }
            
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "error": "Nextflow执行超时（30分钟）",
            "execution_time": time.time() - start_time
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"执行FastP流水线时发生错误: {str(e)}",
            "execution_time": 0
        }


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
        
        return {
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
        
    except Exception as e:
        return {
            "success": False,
            "error": f"解析FastP结果失败: {str(e)}"
        }
