# agent/tools.py - v5.2 Architecture (with all tools implemented)
import os
import time
import json
import threading
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, List
from datetime import datetime


# 移除进度条模块导入
# from .progress_bar import format_progress_update, get_progress_summary

# --- Module-level lock for config file operations ---
_config_lock = threading.Lock()

# --- Environment and Tool Mappings ---
ENVIRONMENT_TOOLS = {
    "qc_env": {
        "fastp": "fastp",
        "fasterq-dump": "sra-tools"
    },
    "align_env": {
        "star": "star"
    },
    "quant_env": {
        "featurecounts": "subread"
    }
}

# --- Tool Dependencies ---
TOOL_DEPENDENCIES = {
    # 环境检查
    "check_environment": [],
    "setup_environment": ["check_environment"],
    
    # 基因组管理
    "search_genome": [],
    "download_genome": ["search_genome"],
    "build_star_index": ["search_genome"],
    
    # 数据管理
    "search_fastq": [],
    "download_fastq": ["search_fastq"],
    "validate_fastq": ["search_fastq"],
    
    # 质量控制
    "run_fastp": ["validate_fastq"],
    
    # 比对分析
    "run_star_align": ["build_star_index", "run_fastp"],
    
    # 定量分析
    "run_featurecounts": ["run_star_align"],
    
    # 结果整理
    "collect_results": ["run_featurecounts"],
    "generate_report": ["run_fastp", "run_featurecounts"]
}

# --- Helper Functions ---
def _get_project_root() -> Path:
    """Returns the absolute path to the project root."""
    return Path(__file__).parent.parent.resolve()

def _get_genomes_config() -> dict:
    """Reads and parses the genomes.json config file."""
    config_path = _get_project_root() / 'config' / 'genomes.json'
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}

# --- Background Threads for Task Execution & Monitoring ---

def _monitor_and_parse_logs(task_id: str, log_file: Path, task_database: dict, db_lock: threading.Lock):
    """
    Tails a log file in a background thread, parses each new line,
    and updates the corresponding step in the task database.
    The thread terminates when the main task is no longer 'running'.
    """
    print(f"[{task_id}] Log monitor thread started for: {log_file}")
    cursor = 0
    while True:
        with db_lock:
            task_status = task_database.get(task_id, {}).get('status')
        
        if task_status not in ['starting', 'running']:
            print(f"[{task_id}] Main task status is '{task_status}'. Stopping log monitor.")
            break

        try:
            if log_file.exists():
                with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                    f.seek(cursor)
                    for line in f:
                        # 简化的日志处理，不再使用复杂的解析
                        cleaned_line = line.strip()
                        with db_lock:
                            if task_id not in task_database: continue
                            # 直接添加到第一个步骤的日志中
                            if task_database[task_id]['details']['steps']:
                                task_database[task_id]['details']['steps'][0]['logs'].append(cleaned_line)
                    cursor = f.tell()
        except Exception as e:
            print(f"[{task_id}] Error reading log file: {e}")

        time.sleep(2) # Poll for new log entries every 2 seconds
    
    print(f"[{task_id}] Log monitor thread finished.")


# --- New Tools for v5.2 Architecture ---

def get_task_status(task_id: str, task_database: dict, db_lock: threading.Lock) -> dict:
    """
    Queries the status of a specific task ID. This is now a simple, thread-safe
    read from the shared task database, as background threads handle all updates.
    """
    print(f"Tool 'get_task_status' called for task: {task_id}")
    with db_lock:
        task_info = task_database.get(task_id)
        if not task_info:
            return {"status": "error", "message": f"Task '{task_id}' not found."}
        
        return {
            "status": "success",
            "task_id": task_id,
            "task_info": json.loads(json.dumps(task_info))  # 确保返回的是可序列化的JSON
        }

def list_available_genomes() -> dict:
    """Reads and returns the content of config/genomes.json."""
    print("Tool 'list_available_genomes' called.")
    genomes_data = _get_genomes_config()
    if not genomes_data:
        return {"error": "Could not read or find 'config/genomes.json'."}
    return {"status": "success", "genomes": genomes_data}


def list_files(path: str, recursive: bool = False) -> dict:
    """
    Lists files and directories in a given path relative to the project root.
    """
    print(f"Tool 'list_files' called for path: {path}, recursive: {recursive}")
    project_root = _get_project_root()
    target_path = (project_root / path).resolve()

    # Security check to prevent directory traversal
    if project_root not in target_path.parents and target_path != project_root:
        return {"status": "error", "message": "Access denied: Path is outside of the project directory."}

    if not target_path.exists():
        return {"status": "error", "message": f"Path does not exist: {target_path}"}
    
    results = []
    try:
        if recursive:
            for root, dirs, files in os.walk(target_path):
                rel_root = Path(root).relative_to(target_path)
                for name in dirs:
                    results.append(str(rel_root / name) + "/")
                for name in files:
                    results.append(str(rel_root / name))
        else:
            for entry in os.listdir(target_path):
                if (target_path / entry).is_dir():
                    results.append(entry + "/")
                else:
                    results.append(entry)
        return {"status": "success", "path": path, "contents": sorted(results)}
    except Exception as e:
        return {"status": "error", "message": str(e)}

def add_genome_to_config(genome_entry_json_str: str) -> dict:
    """
    Adds or updates a complete genome entry in config/genomes.json from a single JSON string.
    The LLM is responsible for generating the full entry, including the top-level key.
    Example input: '{"xenLae2": {"species": "xenopus_laevis", ...}}'
    """
    print("Tool 'add_genome_to_config' called.")

    # 1. 解析由LLM生成的、包含顶级键的JSON字符串
    try:
        new_entry_dict = json.loads(genome_entry_json_str)
        if not isinstance(new_entry_dict, dict) or len(new_entry_dict) != 1:
            raise ValueError("JSON string must represent a dictionary with a single top-level key (the genome version).")
    except (json.JSONDecodeError, ValueError) as e:
        return {"status": "error", "message": f"Invalid JSON string provided: {e}"}

    # 2. 更新配置文件
    config_path = _get_project_root() / 'config' / 'genomes.json'
    with _config_lock:
        try:
            genomes_config = _get_genomes_config()
            
            # 将新的字典条目合并到现有配置中
            genomes_config.update(new_entry_dict)
            
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, 'w') as f:
                # 保持与原文件一致的缩进格式
                json.dump(genomes_config, f, indent=2)
                
            version_key = list(new_entry_dict.keys())[0]
            return {"status": "success", "message": f"Genome '{version_key}' was successfully added/updated in config/genomes.json."}
        except Exception as e:
            return {"status": "error", "message": f"Failed to write to config file: {e}"}

def download_genome_files(genome_name: str, task_database: dict, db_lock: threading.Lock, task_id_counter: int) -> dict:
    """
    A dedicated tool to download genome files and prepare the STAR index.
    It creates a specific plan and reuses the main execution logic.
    NOTE: This tool requires server-injected dependencies.
    """
    print(f"Tool 'download_genome_files' called for: {genome_name}")

    genomes_config = _get_genomes_config()
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {"status": "error", "message": f"Genome '{genome_name}' not found in config/genomes.json. Please add it first using 'add_genome_to_config'."}

    plan = {
        "srr_ids": [],
        "genome_name": genome_name,
        "tools_to_execute": ['download_genome', 'build_star_index'],
        "messages": [
            f"ℹ️ **计划**: 将下载 '{genome_name}' 的基因组源文件。",
            f"ℹ️ **计划**: 将为 '{genome_name}' 新建 STAR 索引。"
        ],
        "is_executable": True,
        "genome_info": genome_info
    }
    
    description = f"Download and prepare genome: {genome_name}"

    return execute_planned_task(plan, description, task_database, db_lock, task_id_counter)

# --- Fallback Tool ---

def unsupported_request(user_request: str) -> dict:
    """
    Handles requests that are not supported by the current tool set.
    """
    return {
        "status": "error",
        "message": f"不支持此请求: {user_request}",
        "suggestions": [
            "请使用 'plan_analysis_task' 来制定分析计划",
            "请使用 'execute_planned_task' 来执行分析任务",
            "请使用 'list_available_genomes' 来查看可用基因组"
        ]
    }

# --- New Search Tools ---

def search_genome_tool(genome_name: str = None, **kwargs) -> dict:
    """
    搜索和列出基因组信息 - 可查询特定基因组或列出所有可用基因组
    """
    print(f"Tool 'search_genome_tool' called for genome '{genome_name}'.")
    project_root = _get_project_root()
    genomes_config = _get_genomes_config()
    
    # 如果没有指定基因组名称，列出所有可用基因组
    if not genome_name:
        available_genomes = []
        for name, info in genomes_config.items():
            fasta_path = project_root / info['fasta']
            gtf_path = project_root / info['gtf']
            star_index_path = fasta_path.parent / 'star_index'
            
            files_exist = fasta_path.is_file() and gtf_path.is_file()
            star_index_exists = star_index_path.is_dir() and (star_index_path / 'SA').is_file()
            
            available_genomes.append({
                "name": name,
                "species": info['species'],
                "version": info['version'],
                "exists_in_filesystem": files_exist,
                "star_index_exists": star_index_exists,
                "fasta_path": str(fasta_path),
                "gtf_path": str(gtf_path),
                "star_index_path": str(star_index_path)
            })
        
        return {
            "status": "success",
            "message": f"找到 {len(available_genomes)} 个可用基因组",
            "available_genomes": available_genomes,
            "total_count": len(available_genomes)
        }
    
    # 检查配置中是否存在
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "not_found",
            "message": f"基因组 '{genome_name}' 在配置中未找到",
            "genome_name": genome_name,
            "exists_in_config": False,
            "exists_in_filesystem": False
        }
    
    # 检查文件系统中是否存在
    fasta_path = project_root / genome_info['fasta']
    gtf_path = project_root / genome_info['gtf']
    star_index_path = fasta_path.parent / 'star_index'
    
    files_exist = fasta_path.is_file() and gtf_path.is_file()
    star_index_exists = star_index_path.is_dir() and (star_index_path / 'SA').is_file()
    
    return {
        "status": "success",
        "genome_name": genome_name,
        "exists_in_config": True,
        "exists_in_filesystem": files_exist,
        "star_index_exists": star_index_exists,
        "fasta_path": str(fasta_path),
        "gtf_path": str(gtf_path),
        "star_index_path": str(star_index_path),
        "genome_info": genome_info
    }

def search_fastq_tool(srr_id: str, **kwargs) -> dict:
    """
    搜索和验证FASTQ文件 - 检查文件存在性并验证完整性
    """
    print(f"Tool 'search_fastq_tool' called for SRR '{srr_id}'.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    
    # 检查 FASTQ 文件
    fastq_dir = data_dir / 'fastq'
    r1_path = fastq_dir / f"{srr_id}_1.fastq.gz"
    r2_path = fastq_dir / f"{srr_id}_2.fastq.gz"
    
    r1_exists = r1_path.is_file()
    r2_exists = r2_path.is_file()
    
    # 验证文件完整性
    validation_results = {}
    if r1_exists:
        try:
            # 检查文件大小
            r1_size = r1_path.stat().st_size
            validation_results["r1_size"] = r1_size
            validation_results["r1_valid"] = r1_size > 0
        except Exception as e:
            validation_results["r1_valid"] = False
            validation_results["r1_error"] = str(e)
    
    if r2_exists:
        try:
            # 检查文件大小
            r2_size = r2_path.stat().st_size
            validation_results["r2_size"] = r2_size
            validation_results["r2_valid"] = r2_size > 0
        except Exception as e:
            validation_results["r2_valid"] = False
            validation_results["r2_error"] = str(e)
    
    # 检查是否为配对端数据
    is_paired = r1_exists and r2_exists
    is_single = r1_exists and not r2_exists
    
    # 综合验证结果
    overall_valid = (r1_exists and validation_results.get("r1_valid", False)) and \
                   (not is_paired or (r2_exists and validation_results.get("r2_valid", False)))
    
    return {
        "status": "success",
        "srr_id": srr_id,
        "r1_exists": r1_exists,
        "r2_exists": r2_exists,
        "r1_path": str(r1_path),
        "r2_path": str(r2_path),
        "is_paired": is_paired,
        "is_single": is_single,
        "validation_results": validation_results,
        "overall_valid": overall_valid,
        "message": f"FASTQ文件验证完成 - {'有效' if overall_valid else '存在问题'}"
    }

# --- Environment Management Tools ---

def check_environment_tool(**kwargs) -> dict:
    """
    检查当前环境状态 - 修复版本
    """
    print("Tool 'check_environment_tool' called.")
    
    # 检查关键工具 - 在conda环境中查找
    tools_status = {}
    
    # 定义工具和对应的conda环境
    key_tools = {
        "fastp": {"env": "qc_env", "description": "质量控制工具"},
        "STAR": {"env": "align_env", "description": "比对工具"}, 
        "featureCounts": {"env": "quant_env", "description": "定量工具"},
        "samtools": {"env": "align_env", "description": "BAM处理工具"},
        "ascp": {"env": None, "description": "Aspera传输工具"}
    }
    
    for tool, info in key_tools.items():
        try:
            if info["env"] is None:
                # 对于nextflow，检查系统PATH
                result = subprocess.run(["which", tool], 
                                      capture_output=True, 
                                      text=True)
                available = result.returncode == 0
                path = result.stdout.strip() if result.returncode == 0 else None
            else:
                # 对于其他工具，检查conda环境中的路径
                env_path = f"/opt/conda/envs/{info['env']}/bin/{tool}"
                available = os.path.exists(env_path)
                path = env_path if available else None
                
                # 如果工具不存在，尝试使用conda run检查
                if not available:
                    try:
                        result = subprocess.run([
                            "bash", "-c", 
                            f"source /opt/conda/bin/activate {info['env']} && which {tool}"
                        ], capture_output=True, text=True)
                        available = result.returncode == 0
                        path = result.stdout.strip() if result.returncode == 0 else None
                    except Exception:
                        pass
            
            tools_status[tool] = {
                "available": available,
                "description": info["description"],
                "path": path,
                "environment": info["env"]
            }
            print(f"工具 {tool}: {'可用' if available else '不可用'} (环境: {info['env']})")
        except Exception as e:
            tools_status[tool] = {
                "available": False,
                "description": info["description"],
                "error": str(e),
                "environment": info["env"]
            }
            print(f"检查工具 {tool} 时出错: {e}")
    
    # 检查conda环境
    conda_envs = ["qc_env", "align_env", "quant_env"]
    env_status = {}
    
    for env in conda_envs:
        try:
            # 检查环境目录是否存在
            env_path = f"/opt/conda/envs/{env}"
            env_status[env] = {
                "exists": os.path.exists(env_path),
                "path": env_path
            }
            print(f"环境 {env}: {'存在' if os.path.exists(env_path) else '不存在'}")
        except Exception as e:
            env_status[env] = {
                "exists": False,
                "error": str(e)
            }
            print(f"检查环境 {env} 时出错: {e}")
    
    # 检查关键目录
    directories = {
        "data_dir": "/app/data",
        "results_dir": "/app/data/results", 
        "fastq_dir": "/app/data/fastq",
        "genome_dir": "/app/genomes"
    }
    
    dir_status = {}
    for name, path in directories.items():
        dir_status[name] = {
            "exists": os.path.exists(path),
            "path": path,
            "writable": os.access(path, os.W_OK) if os.path.exists(path) else False
        }
    
    return {
        "status": "success",
        "message": "环境检查完成",
        "tools_status": tools_status,
        "conda_environments": env_status,
        "directories": dir_status,
        "summary": {
            "total_tools": len(tools_status),
            "available_tools": sum(1 for t in tools_status.values() if t["available"]),
            "total_envs": len(conda_envs),
            "existing_envs": sum(1 for e in env_status.values() if e["exists"])
        }
    }

def setup_environment_tool(environment_type: str = "conda", **kwargs) -> dict:
    """
    设置分析环境 - 修复版本
    """
    print(f"Tool 'setup_environment_tool' called for environment type '{environment_type}'.")
    
    if environment_type != "conda":
        return {
            "status": "error",
            "message": f"不支持的环境类型: {environment_type}，目前只支持 conda"
        }
    
    # 检查conda环境 - 修复版本
    envs_to_check = ["qc_env", "align_env", "quant_env"]
    env_status = {}
    
    for env in envs_to_check:
        try:
            # 检查环境目录是否存在
            env_path = f"/opt/conda/envs/{env}"
            if os.path.exists(env_path):
                # 检查环境中的关键工具
                tools_to_check = {
                    "qc_env": ["fastp"],
                    "align_env": ["STAR", "samtools"],
                    "quant_env": ["featureCounts"]
                }
                
                available_tools = []
                missing_tools = []
                
                for tool in tools_to_check.get(env, []):
                    # 首先检查工具文件是否存在
                    tool_path = f"{env_path}/bin/{tool}"
                    if os.path.exists(tool_path):
                        available_tools.append(tool)
                    else:
                        # 如果文件不存在，尝试使用conda run检查
                        try:
                            result = subprocess.run([
                                "bash", "-c", 
                                f"source /opt/conda/bin/activate {env} && which {tool}"
                            ], capture_output=True, text=True)
                            if result.returncode == 0:
                                available_tools.append(tool)
                            else:
                                missing_tools.append(tool)
                        except Exception:
                            missing_tools.append(tool)
                
                env_status[env] = {
                    "exists": True,
                    "path": env_path,
                    "available_tools": available_tools,
                    "missing_tools": missing_tools,
                    "ready": len(missing_tools) == 0
                }
            else:
                env_status[env] = {
                    "exists": False,
                    "path": env_path,
                    "available_tools": [],
                    "missing_tools": [],
                    "ready": False
                }
        except Exception as e:
            env_status[env] = {
                "exists": False,
                "error": str(e),
                "ready": False
            }
    
    # 创建必要的目录
    directories = [
        "/app/data/results",
        "/app/data/fastq", 
        "/app/data/results/fastp",
        "/app/data/results/star",
        "/app/data/results/featurecounts"
    ]
    
    created_dirs = []
    failed_dirs = []
    
    for dir_path in directories:
        try:
            os.makedirs(dir_path, exist_ok=True)
            created_dirs.append(dir_path)
        except Exception as e:
            failed_dirs.append(f"{dir_path} (错误: {e})")
    
    # 检查conda命令是否可用
    conda_available = False
    conda_path = None
    
    try:
        result = subprocess.run(["which", "conda"], 
                              capture_output=True, 
                              text=True)
        if result.returncode == 0:
            conda_available = True
            conda_path = result.stdout.strip()
    except Exception:
        pass
    
    return {
        "status": "success" if conda_available else "warning",
        "message": "环境设置检查完成",
        "conda_available": conda_available,
        "conda_path": conda_path,
        "environments": env_status,
        "directories": {
            "created": created_dirs,
            "failed": failed_dirs
        },
        "summary": {
            "total_envs": len(envs_to_check),
            "ready_envs": sum(1 for e in env_status.values() if e.get("ready", False)),
            "total_dirs": len(directories),
            "created_dirs": len(created_dirs)
        }
    }


def _parse_fastp_json(json_file_path: str) -> dict:
    """
    解析 fastp JSON 报告文件
    """
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
        
        summary = data.get('summary', {})
        before_filtering = summary.get('before_filtering', {})
        after_filtering = summary.get('after_filtering', {})
        
        return {
            "total_reads": before_filtering.get('total_reads', 0),
            "total_bases": before_filtering.get('total_bases', 0),
            "q20_rate": before_filtering.get('q20_rate', 0),
            "q30_rate": before_filtering.get('q30_rate', 0),
            "clean_reads": after_filtering.get('total_reads', 0),
            "clean_bases": after_filtering.get('total_bases', 0),
            "clean_q20_rate": after_filtering.get('q20_rate', 0),
            "clean_q30_rate": after_filtering.get('q30_rate', 0),
            "filtered_reads": before_filtering.get('total_reads', 0) - after_filtering.get('total_reads', 0),
            "filter_rate": (before_filtering.get('total_reads', 0) - after_filtering.get('total_reads', 0)) / before_filtering.get('total_reads', 1) * 100
        }
    except Exception as e:
        return {"error": f"解析fastp JSON失败: {e}"}

def _parse_star_log(log_file_path: str) -> dict:
    """
    解析 STAR log 文件
    """
    try:
        with open(log_file_path, 'r') as f:
            lines = f.readlines()
        
        result = {}
        for line in lines:
            line = line.strip()
            if "Number of input reads" in line:
                result["input_reads"] = int(line.split()[-1])
            elif "Uniquely mapped reads number" in line:
                result["unique_mapped"] = int(line.split()[-1])
            elif "Number of reads mapped to multiple loci" in line:
                result["multi_mapped"] = int(line.split()[-1])
            elif "Number of reads unmapped" in line:
                result["unmapped"] = int(line.split()[-1])
            elif "Uniquely mapped reads %" in line:
                result["unique_mapping_rate"] = float(line.split()[-1].replace('%', ''))
            elif "Number of reads mapped to too many loci" in line:
                result["too_many_loci"] = int(line.split()[-1])
            elif "Number of reads unmapped: too many mismatches" in line:
                result["too_many_mismatches"] = int(line.split()[-1])
            elif "Number of reads unmapped: too short" in line:
                result["too_short"] = int(line.split()[-1])
            elif "Number of reads unmapped: other" in line:
                result["unmapped_other"] = int(line.split()[-1])
        
        if "input_reads" in result and result["input_reads"] > 0:
            result["total_mapping_rate"] = (result.get("unique_mapped", 0) + result.get("multi_mapped", 0)) / result["input_reads"] * 100
        
        return result
    except Exception as e:
        return {"error": f"解析STAR log失败: {e}"}

def _parse_featurecounts_summary(summary_file_path: str) -> dict:
    """
    解析 featureCounts summary 文件
    """
    try:
        with open(summary_file_path, 'r') as f:
            lines = f.readlines()
        
        result = {}
        for line in lines:
            line = line.strip()
            if line and '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    status = parts[0]
                    count = int(parts[1])
                    result[status] = count
        
        # 计算总reads和assigned比例
        total_reads = sum(result.values())
        assigned_reads = result.get("Assigned", 0)
        
        if total_reads > 0:
            result["assigned_rate"] = assigned_reads / total_reads * 100
            result["total_reads"] = total_reads
        
        return result
    except Exception as e:
        return {"error": f"解析featureCounts summary失败: {e}"}

def collect_results_tool(srr_ids: List[str], **kwargs) -> dict:
    """
    收集分析结果并解析关键指标
    """
    print(f"Tool 'collect_results_tool' called for SRRs {srr_ids}.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    
    # 收集各种结果文件
    results = {
        "fastp_results": [],
        "star_results": [],
        "featurecounts_results": {},
        "analysis_summary": {}
    }
    
    # 解析每个样本的结果
    for srr_id in srr_ids:
        sample_summary = {"srr_id": srr_id}
        
        # fastp 结果解析
        fastp_dir = data_dir / 'results' / 'fastp' / srr_id
        json_file = fastp_dir / f"{srr_id}.json"
        if json_file.exists():
            fastp_stats = _parse_fastp_json(str(json_file))
            sample_summary["fastp_stats"] = fastp_stats
            results["fastp_results"].append({
                "srr_id": srr_id,
                "html_report": str(fastp_dir / f"{srr_id}.html"),
                "json_report": str(json_file),
                "stats": fastp_stats
            })
        
        # STAR 结果解析
        bam_dir = data_dir / 'results' / 'bam' / srr_id
        log_file = bam_dir / f"{srr_id}.Log.final.out"
        if log_file.exists():
            star_stats = _parse_star_log(str(log_file))
            sample_summary["star_stats"] = star_stats
            results["star_results"].append({
                "srr_id": srr_id,
                "bam_file": str(bam_dir / f"{srr_id}.Aligned.sortedByCoord.out.bam"),
                "log_file": str(log_file),
                "stats": star_stats
            })
        
        results["analysis_summary"][srr_id] = sample_summary
    
    # featureCounts 结果解析
    featurecounts_dir = data_dir / 'results' / 'featurecounts'
    summary_file = featurecounts_dir / "counts.txt.summary"
    if summary_file.exists():
        featurecounts_stats = _parse_featurecounts_summary(str(summary_file))
        results["featurecounts_results"] = {
            "counts_file": str(featurecounts_dir / "counts.txt"),
            "summary_file": str(summary_file),
            "stats": featurecounts_stats
        }
    
    return {
        "status": "success",
        "message": f"结果收集和解析完成",
        "srr_ids": srr_ids,
        "results": results
    }

def generate_report_tool(srr_ids: List[str], **kwargs) -> dict:
    """
    生成分析报告 - 包含详细的上游分析总结
    """
    print(f"Tool 'generate_report_tool' called for SRR IDs: {srr_ids}")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    results_dir = data_dir / 'results'
    
    try:
        # 收集并解析结果数据
        results_data = collect_results_tool(srr_ids)
        if results_data["status"] != "success":
            return results_data
        
        results = results_data["results"]
        
        # 创建报告目录
        report_dir = results_dir / 'reports'
        report_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成HTML报告
        report_content = f"""
        <html>
        <head>
            <title>RNA-seq 分析报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
                .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .success {{ color: green; }}
                .error {{ color: red; }}
                .warning {{ color: orange; }}
                table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .metric {{ font-weight: bold; color: #333; }}
                .value {{ font-family: monospace; }}
                .summary-box {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin: 10px 0; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>RNA-seq 分析报告</h1>
                <p><strong>分析时间:</strong> {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>分析样本:</strong> {', '.join(srr_ids)}</p>
            </div>
            
            <div class="section">
                <h2>📊 上游分析总结</h2>
        """
        
        # 添加每个样本的详细统计
        for srr_id in srr_ids:
            sample_summary = results["analysis_summary"].get(srr_id, {})
            report_content += f"""
                <div class="summary-box">
                    <h3>样本: {srr_id}</h3>
            """
            
            # fastp 统计
            if "fastp_stats" in sample_summary:
                fastp_stats = sample_summary["fastp_stats"]
                if "error" not in fastp_stats:
                    report_content += f"""
                        <h4>🔍 质量控制 (fastp)</h4>
                        <table>
                            <tr><th>指标</th><th>原始数据</th><th>过滤后</th></tr>
                            <tr><td>总reads数</td><td class="value">{fastp_stats.get('total_reads', 0):,}</td><td class="value">{fastp_stats.get('clean_reads', 0):,}</td></tr>
                            <tr><td>总碱基数</td><td class="value">{fastp_stats.get('total_bases', 0):,}</td><td class="value">{fastp_stats.get('clean_bases', 0):,}</td></tr>
                            <tr><td>Q20比例</td><td class="value">{fastp_stats.get('q20_rate', 0):.2f}%</td><td class="value">{fastp_stats.get('clean_q20_rate', 0):.2f}%</td></tr>
                            <tr><td>Q30比例</td><td class="value">{fastp_stats.get('q30_rate', 0):.2f}%</td><td class="value">{fastp_stats.get('clean_q30_rate', 0):.2f}%</td></tr>
                            <tr><td>过滤reads数</td><td colspan="2" class="value">{fastp_stats.get('filtered_reads', 0):,} ({fastp_stats.get('filter_rate', 0):.2f}%)</td></tr>
                        </table>
                    """
                else:
                    report_content += f"<p class='error'>fastp统计解析失败: {fastp_stats['error']}</p>"
            
            # STAR 统计
            if "star_stats" in sample_summary:
                star_stats = sample_summary["star_stats"]
                if "error" not in star_stats:
                    report_content += f"""
                        <h4>🎯 序列比对 (STAR)</h4>
                        <table>
                            <tr><th>指标</th><th>数值</th></tr>
                            <tr><td>输入reads数</td><td class="value">{star_stats.get('input_reads', 0):,}</td></tr>
                            <tr><td>唯一比对reads</td><td class="value">{star_stats.get('unique_mapped', 0):,} ({star_stats.get('unique_mapping_rate', 0):.2f}%)</td></tr>
                            <tr><td>多重比对reads</td><td class="value">{star_stats.get('multi_mapped', 0):,}</td></tr>
                            <tr><td>总比对率</td><td class="value">{star_stats.get('total_mapping_rate', 0):.2f}%</td></tr>
                            <tr><td>未比对reads</td><td class="value">{star_stats.get('unmapped', 0):,}</td></tr>
                        </table>
                    """
                else:
                    report_content += f"<p class='error'>STAR统计解析失败: {star_stats['error']}</p>"
            
            report_content += "</div>"
        
        # featureCounts 统计
        if "featurecounts_results" in results and "stats" in results["featurecounts_results"]:
            fc_stats = results["featurecounts_results"]["stats"]
            if "error" not in fc_stats:
                report_content += f"""
                    <div class="summary-box">
                        <h3>📈 基因定量 (featureCounts)</h3>
                        <table>
                            <tr><th>状态</th><th>reads数</th><th>比例</th></tr>
                """
                
                total_reads = fc_stats.get("total_reads", 0)
                for status, count in fc_stats.items():
                    if status not in ["total_reads", "assigned_rate"]:
                        percentage = (count / total_reads * 100) if total_reads > 0 else 0
                        report_content += f"<tr><td>{status}</td><td class='value'>{count:,}</td><td class='value'>{percentage:.2f}%</td></tr>"
                
                report_content += f"""
                            <tr><td><strong>总reads</strong></td><td class='value'><strong>{total_reads:,}</strong></td><td class='value'><strong>100%</strong></td></tr>
                        </table>
                        <p><strong>基因分配率:</strong> <span class='value'>{fc_stats.get('assigned_rate', 0):.2f}%</span></p>
                    </div>
                """
            else:
                report_content += f"<p class='error'>featureCounts统计解析失败: {fc_stats['error']}</p>"
        
        report_content += """
            </div>
            
            <div class="section">
                <h2>✅ 分析步骤完成情况</h2>
                <ul>
                    <li class="success">✓ 数据下载和验证</li>
                    <li class="success">✓ 质量控制 (fastp)</li>
                    <li class="success">✓ 序列比对 (STAR)</li>
                    <li class="success">✓ 基因定量 (featureCounts)</li>
                </ul>
            </div>
            
            <div class="section">
                <h2>📁 结果文件位置</h2>
                <ul>
                    <li><strong>质量控制结果:</strong> <code>data/results/fastp/</code></li>
                    <li><strong>比对结果:</strong> <code>data/results/bam/</code></li>
                    <li><strong>定量结果:</strong> <code>data/results/featurecounts/</code></li>
                    <li><strong>分析报告:</strong> <code>data/results/reports/</code></li>
                </ul>
            </div>
            
            <div class="section">
                <h2>💡 质量评估建议</h2>
                <ul>
                    <li><strong>测序质量:</strong> Q30比例应 > 80%，Q20比例应 > 90%</li>
                    <li><strong>比对质量:</strong> 唯一比对率应 > 70%，总比对率应 > 80%</li>
                    <li><strong>定量质量:</strong> 基因分配率应 > 60%</li>
                </ul>
            </div>
        </body>
        </html>
        """
        
        # 保存报告
        report_file = report_dir / f"report_{'_'.join(srr_ids)}.html"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report_content)
        
        return {
            "status": "success",
            "message": "分析报告生成完成，包含详细的上游分析总结",
            "report_file": str(report_file),
            "srr_ids": srr_ids,
            "results_summary": results["analysis_summary"]
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"生成报告时发生错误: {e}"
        }

# --- React模式专用工具 ---

def react_status_tool(**kwargs) -> dict:
    """
    React模式状态查询工具
    """
    print("Tool 'react_status_tool' called.")
    
    return {
        "status": "success",
        "react_mode": True,
        "current_cycle": kwargs.get("current_cycle", 0),
        "max_cycles": kwargs.get("max_cycles", 10),
        "available_tools": list(ENVIRONMENT_TOOLS.keys()) + [
            "search_genome_tool", "search_fastq_tool", "download_genome_tool",
            "download_fastq_tool", "validate_fastq_tool", "build_star_index_tool",
            "run_fastp_tool", "run_star_align_tool", "run_featurecounts_tool",
            "collect_results_tool", "generate_report_tool"
        ],
        "tool_dependencies": TOOL_DEPENDENCIES
    }

def validate_tool_call_format(tool_name: str, arguments: str) -> dict:
    """
    验证工具调用格式是否正确
    """
    print(f"Tool 'validate_tool_call_format' called for tool '{tool_name}'.")
    
    # 检查是否有禁止的格式
    forbidden_patterns = [
        "REDACTED_SPECIAL_TOKEN",
        "<JSON>",
        "function",
        "```"
    ]
    
    errors = []
    for pattern in forbidden_patterns:
        if pattern in arguments:
            errors.append(f"检测到禁止的格式: {pattern}")
    
    if errors:
        return {
            "status": "error",
            "message": "工具调用格式错误",
            "errors": errors,
            "tool_name": tool_name,
            "arguments": arguments
        }
    
    # 尝试解析JSON
    try:
        import json
        parsed_args = json.loads(arguments)
        return {
            "status": "success",
            "message": "工具调用格式正确",
            "tool_name": tool_name,
            "arguments": parsed_args
        }
    except json.JSONDecodeError as e:
        return {
            "status": "error",
            "message": f"JSON解析错误: {e}",
            "tool_name": tool_name,
            "arguments": arguments
        }


def check_files_exist_tool(file_type: str, **kwargs) -> dict:
    """
    通用文件验证工具 - 检查各种类型的文件是否存在
    """
    print(f"Tool 'check_files_exist_tool' called for file_type '{file_type}'.")
    
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    genomes_config = _get_genomes_config()
    
    results = {
        "file_type": file_type,
        "exists": False,
        "details": {},
        "message": ""
    }
    
    try:
        if file_type == "star_index":
            # 检查STAR索引
            genome_name = kwargs.get('genome_name')
            if not genome_name:
                return {
                    "status": "error",
                    "message": "检查STAR索引需要提供genome_name参数"
                }
            
            genome_info = genomes_config.get(genome_name)
            if not genome_info:
                return {
                    "status": "error",
                    "message": f"基因组 '{genome_name}' 在配置中未找到"
                }
            
            star_index_path = project_root / genome_info['fasta'].replace('.fa', '/star_index')
            exists = star_index_path.is_dir() and (star_index_path / 'SA').is_file()
            
            results.update({
                "exists": exists,
                "file_path": str(star_index_path),
                "details": {
                    "index_dir": str(star_index_path),
                    "has_sa_file": (star_index_path / 'SA').is_file() if star_index_path.is_dir() else False
                },
                "message": f"STAR索引检查完成 - {'存在' if exists else '不存在'}"
            })
            
        elif file_type == "gtf_file":
            # 检查GTF文件
            genome_name = kwargs.get('genome_name')
            if not genome_name:
                return {
                    "status": "error",
                    "message": "检查GTF文件需要提供genome_name参数"
                }
            
            genome_info = genomes_config.get(genome_name)
            if not genome_info:
                return {
                    "status": "error",
                    "message": f"基因组 '{genome_name}' 在配置中未找到"
                }
            
            gtf_path = project_root / genome_info['gtf']
            exists = gtf_path.is_file()
            
            results.update({
                "exists": exists,
                "file_path": str(gtf_path),
                "file_size": gtf_path.stat().st_size if exists else 0,
                "message": f"GTF文件检查完成 - {'存在' if exists else '不存在'}"
            })
            
        elif file_type == "bam_files":
            # 检查BAM文件
            srr_ids = kwargs.get('srr_ids')
            if not srr_ids:
                return {
                    "status": "error",
                    "message": "检查BAM文件需要提供srr_ids参数"
                }
            
            bam_results = {}
            all_exist = True
            
            for srr_id in srr_ids:
                bam_file = work_dir / 'results' / 'bam' / srr_id / f"{srr_id}.bam"
                exists = bam_file.is_file()
                bam_results[srr_id] = {
                    "bam_path": str(bam_file),
                    "exists": exists,
                    "file_size": bam_file.stat().st_size if exists else 0
                }
                if not exists:
                    all_exist = False
            
            results.update({
                "exists": all_exist,
                "srr_ids": srr_ids,
                "details": bam_results,
                "message": f"BAM文件检查完成 - {'全部存在' if all_exist else '部分缺失'}"
            })
            
        elif file_type == "fastq_files":
            # 检查FASTQ文件
            srr_ids = kwargs.get('srr_ids')
            if not srr_ids:
                return {
                    "status": "error",
                    "message": "检查FASTQ文件需要提供srr_ids参数"
                }
            
            fastq_results = {}
            all_exist = True
            
            for srr_id in srr_ids:
                r1_path = work_dir / 'fastq' / f"{srr_id}_1.fastq.gz"
                r2_path = work_dir / 'fastq' / f"{srr_id}_2.fastq.gz"
                r1_exists = r1_path.is_file()
                r2_exists = r2_path.is_file()
                
                fastq_results[srr_id] = {
                    "r1_path": str(r1_path),
                    "r2_path": str(r2_path),
                    "r1_exists": r1_exists,
                    "r2_exists": r2_exists,
                    "is_paired": r1_exists and r2_exists
                }
                
                if not (r1_exists and r2_exists):
                    all_exist = False
            
            results.update({
                "exists": all_exist,
                "srr_ids": srr_ids,
                "details": fastq_results,
                "message": f"FASTQ文件检查完成 - {'全部存在' if all_exist else '部分缺失'}"
            })
            
        elif file_type == "fastp_results":
            # 检查FASTP结果文件
            srr_ids = kwargs.get('srr_ids')
            if not srr_ids:
                return {
                    "status": "error",
                    "message": "检查FASTP结果需要提供srr_ids参数"
                }
            
            fastp_results = {}
            all_exist = True
            
            for srr_id in srr_ids:
                html_file = work_dir / 'results' / 'fastp' / srr_id / f"{srr_id}.html"
                json_file = work_dir / 'results' / 'fastp' / srr_id / f"{srr_id}.json"
                html_exists = html_file.is_file()
                json_exists = json_file.is_file()
                
                fastp_results[srr_id] = {
                    "html_path": str(html_file),
                    "json_path": str(json_file),
                    "html_exists": html_exists,
                    "json_exists": json_exists,
                    "complete": html_exists and json_exists
                }
                
                if not (html_exists and json_exists):
                    all_exist = False
            
            results.update({
                "exists": all_exist,
                "srr_ids": srr_ids,
                "details": fastp_results,
                "message": f"FASTP结果检查完成 - {'全部存在' if all_exist else '部分缺失'}"
            })
            
        elif file_type == "counts_file":
            # 检查featureCounts结果文件
            counts_file = work_dir / 'results' / 'featurecounts' / 'counts.txt'
            exists = counts_file.is_file()
            
            results.update({
                "exists": exists,
                "file_path": str(counts_file),
                "file_size": counts_file.stat().st_size if exists else 0,
                "message": f"featureCounts结果检查完成 - {'存在' if exists else '不存在'}"
            })
            
        else:
            return {
                "status": "error",
                "message": f"不支持的文件类型: {file_type}。支持的类型: star_index, gtf_file, bam_files, fastq_files"
            }
        
        return {
            "status": "success",
            **results
        }
        
    except Exception as e:
        return {
            "status": "error",
            "message": f"检查文件时发生错误: {e}"
        }

def start_analysis_tool(plan: str) -> dict:
    """
    确认分析计划并触发分析状态的转换。
    这个工具是一个信号，告诉服务器从 CONVERSING 切换到 ANALYZING 状态。
    """
    print(f"Tool 'start_analysis_tool' called with plan: {plan}")
    return {
        "status": "success",
        "message": f"分析计划已确认，准备开始执行。计划详情: {plan}"
    }

def run_nextflow_pipeline(
    srr_ids: str = "",
    local_genome_path: str = "",
    local_gtf_path: str = "",
    download_genome_url: str = "",
    download_gtf_url: str = "",
    fastq_dir: str = "",
    fastq_pattern: str = "*_{1,2}.fastq.gz",
    data: str = "./data",
    run_download_srr: bool = False,
    run_download_genome: bool = False,
    run_build_star_index: bool = False,
    run_fastp: bool = False,
    run_star_align: bool = False,
    run_featurecounts: bool = False,
    resume: bool = False,
    star_overhang: int = 100,
    star_threads: int = 4,
    fastp_threads: int = 4,
    featurecounts_threads: int = 4,
    **kwargs
) -> dict:
    """
    运行 Nextflow RNA-seq 分析工作流。
    此函数动态构建并执行 'nextflow run main.nf' 命令。
    """
    print("Tool 'run_nextflow_pipeline' called.")
    project_root = _get_project_root()
    main_nf_path = project_root / 'main.nf'

    if not main_nf_path.exists():
        return {
            "status": "error",
            "message": f"Nextflow 主文件未找到: {main_nf_path}"
        }

    # 构建 Nextflow 命令
    command = ["nextflow", "run", str(main_nf_path), "-with-trace"]

    # --- 添加所有参数 ---
    if srr_ids:
        command.extend(["--srr_ids", srr_ids])
    if local_genome_path:
        command.extend(["--local_genome_path", local_genome_path])
    if local_gtf_path:
        command.extend(["--local_gtf_path", local_gtf_path])
    if download_genome_url:
        command.extend(["--download_genome_url", download_genome_url])
    if download_gtf_url:
        command.extend(["--download_gtf_url", download_gtf_url])
    if fastq_dir:
        command.extend(["--fastq_dir", fastq_dir])
    if fastq_pattern:
        command.extend(["--fastq_pattern", fastq_pattern])
    if data:
        command.extend(["--data", data])

    # --- 添加布尔控制参数 ---
    if run_download_srr:
        command.append("--run_download_srr")
    if run_download_genome:
        command.append("--run_download_genome")
    if run_build_star_index:
        command.append("--run_build_star_index")
    if run_fastp:
        command.append("--run_fastp")
    if run_star_align:
        command.append("--run_star_align")
    if run_featurecounts:
        command.append("--run_featurecounts")
    
    # --- 添加其他参数 ---
    if resume:
        command.append("-resume") # Nextflow resume flag
    command.extend(["--star_overhang", str(star_overhang)])
    command.extend(["--star_threads", str(star_threads)])
    command.extend(["--fastp_threads", str(fastp_threads)])
    command.extend(["--featurecounts_threads", str(featurecounts_threads)])

    print(f"执行 Nextflow 命令: {' '.join(command)}")

    try:
        # 设置环境变量
        env = os.environ.copy()
        # 确保 NXF_HOME 设置在项目目录内，以避免权限问题
        env['NXF_HOME'] = str(project_root / '.nextflow')
        
        # 执行命令
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            env=env,
            cwd=project_root  # 在项目根目录下运行
        )

        if result.returncode == 0:
            return {
                "status": "success",
                "message": "Nextflow 工作流执行成功",
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": " ".join(command),
            }
        else:
            return {
                "status": "error",
                "message": "Nextflow 工作流执行失败",
                "stdout": result.stdout,
                "stderr": result.stderr,
                "command": " ".join(command),
            }
    except FileNotFoundError:
        return {
            "status": "error",
            "message": "Nextflow 命令未找到。请确保 Nextflow 已安装并在 PATH 中。"
        }
    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "message": "Nextflow 工作流执行超时。"
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"运行 Nextflow 工作流时发生错误: {e}"
        }

def _setup_pipeline_logger() -> logging.Logger:
    """
    设置 RNA-seq 流程专用的日志记录器
    """
    project_root = _get_project_root()
    log_file = project_root / 'para.log.txt'
    
    # 创建日志记录器
    logger = logging.getLogger('rna_seq_pipeline')
    logger.setLevel(logging.INFO)
    
    # 避免重复添加处理器
    if not logger.handlers:
        # 创建文件处理器
        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(logging.INFO)
        
        # 创建格式化器
        formatter = logging.Formatter(
            '%(asctime)s - [RNA-seq Pipeline] - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        file_handler.setFormatter(formatter)
        
        # 添加处理器到日志记录器
        logger.addHandler(file_handler)
    
    return logger

def execute_rna_seq_pipeline(
    srr_ids: str = "",
    local_genome_path: str = "",
    local_gtf_path: str = "",
    download_genome_url: str = "",
    download_gtf_url: str = "",
    fastq_dir: str = "",
    fastq_pattern: str = "*_{1,2}.fastq.gz",
    data: str = "./data",
    run_download_srr: bool = False,
    run_download_genome: bool = False,
    run_build_star_index: bool = False,
    run_fastp: bool = False,
    run_star_align: bool = False,
    run_featurecounts: bool = False,
    resume: bool = False,
    star_overhang: int = 100,
    star_threads: int = 4,
    fastp_threads: int = 4,
    featurecounts_threads: int = 4,
    **kwargs
) -> dict:
    """
    在 ANALYZING 状态下，运行完整的 RNA-seq 分析工作流。
    此工具接收所有必要的参数，并一次性启动 Nextflow 流程。
    """
    # 设置日志记录器
    logger = _setup_pipeline_logger()
    
    # 记录函数开始执行和参数信息
    logger.info("=" * 80)
    logger.info("开始执行 RNA-seq 分析流程")
    logger.info(f"执行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("输入参数:")
    logger.info(f"  - srr_ids: {srr_ids}")
    logger.info(f"  - local_genome_path: {local_genome_path}")
    logger.info(f"  - local_gtf_path: {local_gtf_path}")
    logger.info(f"  - download_genome_url: {download_genome_url}")
    logger.info(f"  - download_gtf_url: {download_gtf_url}")
    logger.info(f"  - fastq_dir: {fastq_dir}")
    logger.info(f"  - fastq_pattern: {fastq_pattern}")
    logger.info(f"  - data: {data}")
    logger.info(f"  - run_download_srr: {run_download_srr}")
    logger.info(f"  - run_download_genome: {run_download_genome}")
    logger.info(f"  - run_build_star_index: {run_build_star_index}")
    logger.info(f"  - run_fastp: {run_fastp}")
    logger.info(f"  - run_star_align: {run_star_align}")
    logger.info(f"  - run_featurecounts: {run_featurecounts}")
    logger.info(f"  - resume: {resume}")
    logger.info(f"  - star_overhang: {star_overhang}")
    logger.info(f"  - star_threads: {star_threads}")
    logger.info(f"  - fastp_threads: {fastp_threads}")
    logger.info(f"  - featurecounts_threads: {featurecounts_threads}")
    
    print(f"Tool 'execute_rna_seq_pipeline' called.")
    
    try:
        # 记录开始调用 Nextflow 工作流
        logger.info("开始调用 run_nextflow_pipeline 函数")
        
        # 调用 Nextflow 工作流
        result = run_nextflow_pipeline(
            srr_ids=srr_ids,
            local_genome_path=local_genome_path,
            local_gtf_path=local_gtf_path,
            download_genome_url=download_genome_url,
            download_gtf_url=download_gtf_url,
            fastq_dir=fastq_dir,
            fastq_pattern=fastq_pattern,
            data=data,
            run_download_srr=run_download_srr,
            run_download_genome=run_download_genome,
            run_build_star_index=run_build_star_index,
            run_fastp=run_fastp,
            run_star_align=run_star_align,
            run_featurecounts=run_featurecounts,
            resume=resume,
            star_overhang=star_overhang,
            star_threads=star_threads,
            fastp_threads=fastp_threads,
            featurecounts_threads=featurecounts_threads,
            **kwargs
        )
        
        # 记录执行结果
        if result.get('status') == 'success':
            logger.info("✅ run_nextflow_pipeline 执行成功")
            logger.info(f"执行的命令: {result.get('command', 'N/A')}")
            if result.get('stdout'):
                logger.info("标准输出:")
                for line in result['stdout'].split('\n'):
                    if line.strip():
                        logger.info(f"  {line}")
        else:
            logger.error("❌ run_nextflow_pipeline 执行失败")
            logger.error(f"错误信息: {result.get('message', 'N/A')}")
            logger.error(f"执行的命令: {result.get('command', 'N/A')}")
            if result.get('stderr'):
                logger.error("标准错误输出:")
                for line in result['stderr'].split('\n'):
                    if line.strip():
                        logger.error(f"  {line}")
        
        # 记录函数完成状态
        logger.info(f"RNA-seq 分析流程执行完成，状态: {result.get('status', 'unknown')}")
        logger.info("=" * 80)
        
        return result
        
    except Exception as e:
        # 记录异常信息
        logger.error(f"❌ 执行 RNA-seq 分析流程时发生异常: {str(e)}")
        logger.error(f"异常类型: {type(e).__name__}")
        logger.error("=" * 80)
        
        # 重新抛出异常，保持原有的错误处理逻辑
        raise e