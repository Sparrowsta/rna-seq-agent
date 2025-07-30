# agent/tools.py - v5.2 Architecture (with all tools implemented)
import os
import time
import json
import threading
import subprocess
from pathlib import Path
from typing import Dict, Any, List

# --- New Imports ---


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
                        step_name, cleaned_line = parse_log_line(line)
                        with db_lock:
                            if task_id not in task_database: continue
                            for step in task_database[task_id]['details']['steps']:
                                if step['name'] == step_name:
                                    if step['status'] == 'pending':
                                        step['status'] = 'running'
                                    step['logs'].append(cleaned_line)
                                    break
                    cursor = f.tell()
        except Exception as e:
            print(f"[{task_id}] Error reading log file: {e}")

        time.sleep(2) # Poll for new log entries every 2 seconds
    
    print(f"[{task_id}] Log monitor thread finished.")

def _execute_task(task_id: str, plan: dict, task_database: dict, db_lock: threading.Lock):
    """
    The main background thread that executes a planned Nextflow task.
    It builds the command, starts the process and the log monitor,
    and updates the final status upon completion.
    """
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    log_file = work_dir / 'logs' / f"{task_id}.log"

    try:
        # 1. Build Nextflow command from the execution plan
        command = ["nextflow", "run", str(project_root / "nextflow" / "main.nf"), "-resume"]
        
        # 将工具级计划转换为 Nextflow 步骤
        tools_to_steps = {
            'download_genome': 'genome_download',
            'build_star_index': 'genome_prep',
            'download_fastq': 'download_fastq',
            'run_fastp': 'analysis',
            'run_star_align': 'analysis',
            'run_featurecounts': 'analysis',
            'validate_fastq': 'analysis'
        }
        
        steps_to_execute = []
        for tool in plan.get('tools_to_execute', []):
            if tool in tools_to_steps:
                step = tools_to_steps[tool]
                if step not in steps_to_execute:
                    steps_to_execute.append(step)
        
        for step in steps_to_execute:
            command.extend([f"--run_{step.lower()}", "true"])

        # The main.nf script is smart enough to look up genome info from the config
        # based on the name. We only need to provide the name.
        command.extend(["--genome_name", plan['genome_name']])

        # --- FIX: Re-fetch full genome_info to ensure 'species' and 'version' are present ---
        # This makes _execute_task robust against LLM potentially truncating the plan['genome_info']
        full_genome_config = _get_genomes_config()
        genome_info = full_genome_config.get(plan['genome_name'])
        if not genome_info:
            raise ValueError(f"Genome '{plan['genome_name']}' not found in config/genomes.json during task execution.")

        # Add species and genome_version parameters for proper directory structure
        command.extend(["--species", genome_info['species']])
        command.extend(["--genome_version", genome_info['version']])

        # --- ENHANCED FIX: Check if genome files exist and add download step if needed ---
        # This makes the system robust against LLM truncating the plan['tools_to_execute']
        fasta_path = project_root / genome_info['fasta']
        gtf_path = project_root / genome_info['gtf']
        files_exist = fasta_path.is_file() and gtf_path.is_file()
        
        # Only add genome_download if files don't exist AND it's not already in the plan
        if not files_exist and 'genome_download' not in steps_to_execute:
            print(f"[{task_id}] WARNING: Genome files don't exist but 'genome_download' not in plan. Adding it.")
            steps_to_execute.insert(0, 'genome_download')  # Add at beginning
            command.insert(-2, "--run_genome_download")
            command.insert(-1, "true")
        elif files_exist and 'genome_download' in steps_to_execute:
            print(f"[{task_id}] INFO: Genome files exist, removing 'genome_download' from plan.")
            steps_to_execute.remove('genome_download')  # Remove if files exist

        # Only provide FASTA/GTF paths if we are NOT downloading them,
        # meaning they must already exist.
        if 'genome_download' not in steps_to_execute:
            command.extend(["--fasta", str(project_root / genome_info['fasta'])])
            command.extend(["--gtf", str(project_root / genome_info['gtf'])])

        # Provide SRR IDs for download, or a glob pattern for existing reads.
        if 'download_fastq' in steps_to_execute:
            command.extend(["--srr_ids", ",".join(plan['srr_ids'])])
        elif plan['srr_ids']: # Only add --reads if there are SRR IDs to analyze
            srr_glob = "{" + ",".join(plan['srr_ids']) + "}"
            reads_pattern = str(work_dir / 'fastq' / f"{srr_glob}_*.fastq.gz")
            command.extend(["--reads", reads_pattern])

        with db_lock:
            task_database[task_id]['status'] = 'running'
            task_database[task_id]['details']['message'] = 'Nextflow process starting...'
            task_database[task_id]['details']['command'] = " ".join(command)

        # 2. Start the log monitor thread
        monitor_thread = threading.Thread(
            target=_monitor_and_parse_logs, args=(task_id, log_file, task_database, db_lock), daemon=True
        )
        monitor_thread.start()

        # 3. Execute the command using our unified executor
        #    CRITICAL: Set NXF_HOME to a writable directory inside the container
        #    to prevent Nextflow from trying to write to the root filesystem.
        nf_env = os.environ.copy()
        nxf_home = work_dir / '.nextflow'
        nxf_home.mkdir(parents=True, exist_ok=True)
        # Ensure the directory has proper permissions
        os.chmod(nxf_home, 0o755)
        nf_env['NXF_HOME'] = str(nxf_home.absolute())

        process = execute_with_logging(command, log_file, work_dir, env=nf_env)
        with db_lock:
            task_database[task_id]['details']['pid'] = process.pid

        # 4. Wait for the process to complete
        exit_code = process.wait()

        # 5. Update final status in the database
        with db_lock:
            final_status = "completed" if exit_code == 0 else "failed"
            task_database[task_id]['status'] = final_status
            task_database[task_id]['end_time'] = time.time()
            task_database[task_id]['details']['exit_code'] = exit_code
            task_database[task_id]['details']['message'] = f"Process finished with exit code {exit_code}."
            for step in task_database[task_id]['details']['steps']:
                if step['status'] == 'running':
                    step['status'] = final_status

    except Exception as e:
        print(f"[{task_id}] Critical error in _execute_task thread: {e}")
        with db_lock:
                task_database[task_id]['status'] = 'failed'
                task_database[task_id]['details']['error'] = str(e)

# --- New Tools for v5.2 Architecture ---

def plan_analysis_task(srr_ids: str, genome_name: str) -> dict:
    """
    制定分析计划 - 工具级版本
    返回详细的工具执行计划，而不是 Nextflow 步骤
    """
    print(f"Tool 'plan_analysis_task' called for genome '{genome_name}' and SRRs '{srr_ids}'.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    
    # 解析 SRR IDs
    srr_list = sorted(list(set(srr_ids.replace(",", " ").split())))
    
    plan = {
        "srr_ids": srr_list,
        "genome_name": genome_name,
        "tools_to_execute": [],
        "messages": [],
        "is_executable": False,
        "genome_info": None,
        "environment_requirements": list(ENVIRONMENT_TOOLS.keys()),
        "estimated_time": len(srr_list) * 30,  # 每个样本30分钟估算
        "dependencies": TOOL_DEPENDENCIES
    }

    # --- Step 1: 验证基因组配置 ---
    genomes_config = _get_genomes_config()
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        plan['messages'].append(f"❌ **错误**: 在配置中找不到基因组 '{genome_name}'。请从可用基因组列表中选择。")
        return {"status": "success", "plan": plan}
    
    plan['genome_info'] = genome_info.copy()

    # --- Step 2: 检查基因组文件 ---
    fasta_path = project_root / plan['genome_info']['fasta']
    gtf_path = project_root / plan['genome_info']['gtf']
    star_index_path = fasta_path.parent / 'star_index'
    
    source_files_exist = fasta_path.is_file() and gtf_path.is_file()
    star_index_exists = star_index_path.is_dir() and (star_index_path / 'SA').is_file()
    
    if not source_files_exist:
        plan['tools_to_execute'].append('download_genome')
        plan['messages'].append(f"ℹ️ **计划**: 将下载 '{genome_name}' 的基因组文件。")
    else:
        plan['messages'].append(f"✅ **基因组文件**: '{genome_name}' 的 FASTA 和 GTF 文件已存在。")
    
    if not star_index_exists:
        plan['tools_to_execute'].append('build_star_index')
        plan['messages'].append(f"ℹ️ **计划**: 将为 '{genome_name}' 构建 STAR 索引。")
    else:
        plan['messages'].append(f"✅ **STAR 索引**: 已存在且完整。")

    # --- Step 3: 检查 FASTQ 文件 ---
    missing_srr_files = []
    for srr_id in srr_list:
        r1_path = data_dir / 'fastq' / f"{srr_id}_1.fastq.gz"
        if not r1_path.exists():
            missing_srr_files.append(srr_id)
    
    if missing_srr_files:
        plan['tools_to_execute'].extend(['download_fastq'] * len(missing_srr_files))
        plan['messages'].append(f"ℹ️ **计划**: 将下载 {len(missing_srr_files)} 个缺失的样本文件: {', '.join(missing_srr_files)}。")
    else:
        plan['messages'].append("✅ **样本文件**: 所有需要的 FASTQ 文件均已存在。")

    # --- Step 4: 添加分析工具 ---
    for srr_id in srr_list:
        plan['tools_to_execute'].extend([
            'validate_fastq',
            'run_fastp',
            'run_star_align',
        ])
    
    # 添加定量和结果工具
    plan['tools_to_execute'].extend([
        'run_featurecounts',
        'collect_results',
        'generate_report'
    ])
    
    plan['messages'].append("➡️ **分析步骤**: 将对所有样本执行质量控制、比对和定量分析。")

    # --- Step 5: 验证计划可行性 ---
    plan['is_executable'] = True
    
    # 计算预估时间
    total_time = 0
    for tool in plan['tools_to_execute']:
        if tool == 'download_genome':
            total_time += 10  # 10分钟
        elif tool == 'build_star_index':
            total_time += 20  # 20分钟
        elif tool == 'download_fastq':
            total_time += 5   # 每个样本5分钟
        elif tool == 'run_fastp':
            total_time += 3   # 每个样本3分钟
        elif tool == 'run_star_align':
            total_time += 15  # 每个样本15分钟
        elif tool == 'run_featurecounts':
            total_time += 5   # 5分钟
        else:
            total_time += 1   # 其他工具1分钟
    
    plan['estimated_time'] = total_time

    return {"status": "success", "plan": plan}

def execute_planned_task(plan: dict, description: str, task_database: dict, db_lock: threading.Lock, task_id_counter: int) -> dict:
    """
    Executes a plan that has been confirmed by the user.
    It creates a new task in the database, initializes the UI step structure,
    and launches the background execution thread.
    """
    print(f"Tool 'execute_planned_task' called.")
    if not plan or not plan.get('is_executable', False):
        return {"status": "error", "message": "The provided plan is invalid or not executable."}

    with db_lock:
        task_id = f"task_{task_id_counter}"
        
        # 将工具级计划转换为 Nextflow 步骤
        tools_to_steps = {
            'download_genome': 'genome_download',
            'build_star_index': 'genome_prep',
            'download_fastq': 'download_fastq',
            'run_fastp': 'analysis',
            'run_star_align': 'analysis',
            'run_featurecounts': 'analysis',
            'validate_fastq': 'analysis'
        }
        
        steps_to_execute = []
        for tool in plan.get('tools_to_execute', []):
            if tool in tools_to_steps:
                step = tools_to_steps[tool]
                if step not in steps_to_execute:
                    steps_to_execute.append(step)
        
        steps = []
        if 'genome_download' in steps_to_execute:
            steps.append({"name": "WGET_GENOME_FILES", "status": "pending", "logs": []})
        if 'genome_prep' in steps_to_execute:
            steps.append({"name": "STAR_GENOME_GENERATE", "status": "pending", "logs": []})
        if 'download_fastq' in steps_to_execute:
            steps.append({"name": "FASTERQ_DUMP", "status": "pending", "logs": []})
        if 'analysis' in steps_to_execute:
            steps.append({"name": "FASTP", "status": "pending", "logs": []})
            steps.append({"name": "STAR_ALIGN", "status": "pending", "logs": []})
            steps.append({"name": "FEATURECOUNTS", "status": "pending", "logs": []})

        task_database[task_id] = {
            "type": "pipeline", "status": "starting", "creation_time": time.time(),
            "start_time": time.time(), "end_time": None,
            "details": {
                "description": description, "message": "Task created, execution thread starting.",
                "plan": plan, "steps": steps, "error": None, "pid": None
            }
        }
    
        threading.Thread(target=_execute_task, args=(task_id, plan, task_database, db_lock), daemon=True).start()

        return {
            "status": "success",
            "message": f"Task {task_id} has been successfully launched.",
            "task_id": task_id, "updated_task_id_counter": task_id_counter + 1
        }

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

def add_genome_to_config(genome_name: str, species: str, fasta_url: str, gtf_url: str) -> dict:
    """
    Adds or updates a genome configuration in config/genomes.json.
    This operation is thread-safe.
    """
    print(f"Tool 'add_genome_to_config' called for: {genome_name}")
    
    config_path = _get_project_root() / 'config' / 'genomes.json'
    
    with _config_lock:
        try:
            genomes_config = _get_genomes_config()
            
            new_entry = {
                "species": species,
                "fasta": f"data/genomes/{species}/{genome_name}/genome.fa",
                "gtf": f"data/genomes/{species}/{genome_name}/genes.gtf",
                "fasta_url": fasta_url,
                "gtf_url": gtf_url
            }
            
            genomes_config[genome_name] = new_entry
            
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, 'w') as f:
                json.dump(genomes_config, f, indent=4)
                
            return {"status": "success", "message": f"Genome '{genome_name}' was successfully added/updated."}
        except Exception as e:
            return {"status": "error", "message": f"Failed to update config: {e}"}

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

def search_genome_tool(genome_name: str, **kwargs) -> dict:
    """
    搜索现有基因组文件
    """
    print(f"Tool 'search_genome_tool' called for genome '{genome_name}'.")
    project_root = _get_project_root()
    genomes_config = _get_genomes_config()
    
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
    搜索现有 FASTQ 文件
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
    
    return {
        "status": "success",
        "srr_id": srr_id,
        "r1_exists": r1_exists,
        "r2_exists": r2_exists,
        "r1_path": str(r1_path),
        "r2_path": str(r2_path),
        "is_paired": r1_exists and r2_exists,
        "is_single": r1_exists and not r2_exists
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
        "nextflow": {"env": None, "description": "工作流管理工具"}
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

# --- Download and Validation Tools ---

def download_genome_tool(genome_name: str, **kwargs) -> dict:
    """
    下载基因组文件
    """
    print(f"Tool 'download_genome_tool' called for genome '{genome_name}'.")
    project_root = _get_project_root()
    genomes_config = _get_genomes_config()
    
    # 检查基因组配置
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "error",
            "message": f"基因组 '{genome_name}' 在配置中未找到"
        }
    
    # 构建下载命令
    fasta_url = genome_info['fasta_url']
    gtf_url = genome_info['gtf_url']
    fasta_path = project_root / genome_info['fasta']
    gtf_path = project_root / genome_info['gtf']
    
    # 创建目录
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    gtf_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # 下载 FASTA 文件
        os.system(f"wget --quiet -O {fasta_path}.gz '{fasta_url}'")
        os.system(f"gunzip -f {fasta_path}.gz")
        
        # 下载 GTF 文件
        os.system(f"wget --quiet -O {gtf_path}.gz '{gtf_url}'")
        os.system(f"gunzip -f {gtf_path}.gz")
        
        return {
            "status": "success",
            "message": f"基因组 '{genome_name}' 下载完成",
            "fasta_path": str(fasta_path),
            "gtf_path": str(gtf_path)
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"下载基因组时发生错误: {e}"
        }

def download_fastq_tool(srr_id: str, **kwargs) -> dict:
    """
    下载 FASTQ 文件
    """
    print(f"Tool 'download_fastq_tool' called for SRR '{srr_id}'.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    fastq_dir = data_dir / 'fastq'
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # 切换到 FASTQ 目录
        os.chdir(fastq_dir)
        
        # 使用 fasterq-dump 下载
        result = os.system(f"fasterq-dump {srr_id} --split-files --gzip -O . -p")
        
        if result == 0:
            return {
                "status": "success",
                "message": f"SRR '{srr_id}' 下载完成",
                "srr_id": srr_id,
                "output_dir": str(fastq_dir)
            }
        else:
            return {
                "status": "error",
                "message": f"下载 SRR '{srr_id}' 失败"
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"下载 FASTQ 时发生错误: {e}"
        }

def validate_fastq_tool(srr_id: str, **kwargs) -> dict:
    """
    验证 FASTQ 文件
    """
    print(f"Tool 'validate_fastq_tool' called for SRR '{srr_id}'.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    fastq_dir = data_dir / 'fastq'
    
    r1_path = fastq_dir / f"{srr_id}_1.fastq.gz"
    r2_path = fastq_dir / f"{srr_id}_2.fastq.gz"
    
    # 检查文件是否存在
    r1_exists = r1_path.is_file()
    r2_exists = r2_path.is_file()
    
    if not r1_exists:
        return {
            "status": "error",
            "message": f"FASTQ 文件不存在: {r1_path}"
        }
    
    # 验证文件完整性
    validation_results = {}
    
    try:
        # 检查 R1 文件
        result = os.system(f"gzip -t {r1_path}")
        validation_results["r1_valid"] = result == 0
        
        # 检查 R2 文件（如果存在）
        if r2_exists:
            result = os.system(f"gzip -t {r2_path}")
            validation_results["r2_valid"] = result == 0
        else:
            validation_results["r2_valid"] = None
        
        return {
            "status": "success",
            "srr_id": srr_id,
            "r1_exists": r1_exists,
            "r2_exists": r2_exists,
            "validation_results": validation_results,
            "is_paired": r1_exists and r2_exists
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"验证 FASTQ 文件时发生错误: {e}"
        }

# --- Build and Analysis Tools ---

def build_star_index_tool(genome_name: str, **kwargs) -> dict:
    """
    构建 STAR 索引 - 修复版本
    """
    print(f"Tool 'build_star_index_tool' called for genome '{genome_name}'.")
    project_root = _get_project_root()
    genomes_config = _get_genomes_config()
    
    # 检查基因组配置
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "error",
            "message": f"基因组 '{genome_name}' 在配置中未找到"
        }
    
    fasta_path = project_root / genome_info['fasta']
    gtf_path = project_root / genome_info['gtf']
    star_index_path = fasta_path.parent / 'star_index'
    
    # 检查输入文件
    if not fasta_path.is_file():
        return {
            "status": "error",
            "message": f"FASTA 文件不存在: {fasta_path}"
        }
    
    if not gtf_path.is_file():
        return {
            "status": "error",
            "message": f"GTF 文件不存在: {gtf_path}"
        }
    
    try:
        # 创建索引目录
        star_index_path.mkdir(parents=True, exist_ok=True)
        
        # 使用subprocess和正确的conda环境激活方式
        command = [
            "bash", "-c",
            f"source /opt/conda/bin/activate align_env && "
            f"STAR --runMode genomeGenerate "
            f"--genomeDir {star_index_path} "
            f"--genomeFastaFiles {fasta_path} "
            f"--sjdbGTFfile {gtf_path} "
            f"--runThreadN 4"
        ]
        
        print(f"执行STAR索引构建命令: {' '.join(command)}")
        
        # 设置环境变量
        env = os.environ.copy()
        env['NXF_HOME'] = str(project_root / 'data' / '.nextflow')
        
        # 执行命令
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        
        if result.returncode == 0:
            return {
                "status": "success",
                "message": f"STAR 索引构建完成",
                "star_index_path": str(star_index_path),
                "stdout": result.stdout
            }
        else:
            return {
                "status": "error",
                "message": f"STAR 索引构建失败: {result.stderr}",
                "stdout": result.stdout,
                "stderr": result.stderr
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"构建 STAR 索引时发生错误: {e}"
        }

def run_fastp_tool(srr_id: str, **kwargs) -> dict:
    """
    运行 fastp 质量控制 - 直接调用独立的fastp.nf脚本
    """
    print(f"Tool 'run_fastp_tool' called for SRR '{srr_id}'.")
    
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    
    # 检查输入文件
    fastq_dir = work_dir / 'fastq'
    r1_path = fastq_dir / f"{srr_id}_1.fastq.gz"
    r2_path = fastq_dir / f"{srr_id}_2.fastq.gz"
    
    if not r1_path.is_file():
        return {
            "status": "error",
            "message": f"输入文件不存在: {r1_path}"
        }
    
    try:
        # 直接调用独立的fastp.nf脚本
        command = [
            "nextflow", "run", str(project_root / "nextflow" / "modules" / "fastp.nf"),
            "--reads", str(fastq_dir / f"{srr_id}_*.fastq.gz"),
            "--outdir", str(work_dir / 'results'),
            "-resume"
        ]
        
        # 设置环境变量
        env = os.environ.copy()
        env['NXF_HOME'] = str(work_dir / '.nextflow')
        
        print(f"执行Nextflow FASTP命令: {' '.join(command)}")
        
        # 执行命令
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        
        if result.returncode == 0:
            # 检查输出文件
            output_dir = work_dir / 'results' / 'fastp' / srr_id
            html_file = output_dir / f"{srr_id}.html"
            json_file = output_dir / f"{srr_id}.json"
            
            return {
                "status": "success",
                "message": f"fastp 质量控制完成",
                "srr_id": srr_id,
                "output_dir": str(output_dir),
                "html_report": str(html_file) if html_file.exists() else None,
                "json_report": str(json_file) if json_file.exists() else None,
                "stdout": result.stdout
            }
        else:
            return {
                "status": "error",
                "message": f"fastp 质量控制失败: ",
                "stdout": result.stdout,
                "stderr": result.stderr
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"运行 fastp 时发生错误: {e}"
        }

def run_star_align_tool(srr_id: str, genome_name: str, **kwargs) -> dict:
    """
    运行 STAR 比对 - 直接调用独立的star.nf脚本
    """
    print(f"Tool 'run_star_align_tool' called for SRR '{srr_id}' with genome '{genome_name}'.")
    
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    genomes_config = _get_genomes_config()
    
    # 检查基因组配置
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "error",
            "message": f"基因组 '{genome_name}' 在配置中未找到"
        }
    
    # 获取STAR索引路径（由LLM决定是否已存在）
    star_index_path = project_root / genome_info['fasta'].replace('.fa', '/star_index')
    
    try:
        # 直接调用独立的star.nf脚本
        command = [
            "nextflow", "run", str(project_root / "nextflow" / "modules" / "star.nf"),
            "--reads", str(work_dir / 'fastq' / f"{srr_id}_*.fastq.gz"),
            "--star_index", str(star_index_path),
            "--outdir", str(work_dir / 'results'),
            "-resume"
        ]
        
        # 设置环境变量
        env = os.environ.copy()
        env['NXF_HOME'] = str(work_dir / '.nextflow')
        
        print(f"执行Nextflow STAR命令: {' '.join(command)}")
        
        # 执行命令
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        
        if result.returncode == 0:
            # 检查输出文件
            output_dir = work_dir / 'results' / 'bam' / srr_id
            bam_file = output_dir / f"{srr_id}.bam"
            log_file = output_dir / f"{srr_id}.log"
            
            return {
                "status": "success",
                "message": f"STAR 比对完成",
                "srr_id": srr_id,
                "output_dir": str(output_dir),
                "bam_file": str(bam_file) if bam_file.exists() else None,
                "log_file": str(log_file) if log_file.exists() else None,
                "stdout": result.stdout
            }
        else:
            return {
                "status": "error",
                "message": f"STAR 比对失败: ",
                "stdout": result.stdout,
                "stderr": result.stderr
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"运行 STAR 比准时发生错误: {e}"
        }

# --- Quantification and Results Tools ---

def run_featurecounts_tool(srr_ids: List[str], genome_name: str, **kwargs) -> dict:
    """
    运行 featureCounts 定量分析 - 直接调用独立的featurecounts.nf脚本
    """
    print(f"Tool 'run_featurecounts_tool' called for SRRs {srr_ids} with genome '{genome_name}'.")
    
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    genomes_config = _get_genomes_config()
    
    # 检查基因组配置
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "error",
            "message": f"基因组 '{genome_name}' 在配置中未找到"
        }
    
    # 获取GTF文件路径（由LLM决定是否已存在）
    gtf_path = project_root / genome_info['gtf']
    
    # 获取BAM文件路径（由LLM决定是否已存在）
    bam_files = []
    for srr_id in srr_ids:
        bam_file = work_dir / 'results' / 'bam' / srr_id / f"{srr_id}.bam"
        bam_files.append(str(bam_file))
    
    try:
        # 直接调用独立的featurecounts.nf脚本
        command = [
            "nextflow", "run", str(project_root / "nextflow" / "modules" / "featurecounts.nf"),
            "--bams", ",".join(bam_files),
            "--gtf", str(gtf_path),
            "--outdir", str(work_dir / 'results'),
            "-resume"
        ]
        
        # 设置环境变量
        env = os.environ.copy()
        env['NXF_HOME'] = str(work_dir / '.nextflow')
        
        print(f"执行Nextflow FEATURECOUNTS命令: {' '.join(command)}")
        
        # 执行命令
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        
        if result.returncode == 0:
            return {
                "status": "success",
                "message": f"featureCounts 定量分析完成",
                "srr_ids": srr_ids,
                "output_dir": str(work_dir / 'results' / 'featurecounts')
            }
        else:
            return {
                "status": "error",
                "message": f"featureCounts 定量分析失败: ",
                "stdout": result.stdout,
                "stderr": result.stderr
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"运行 featureCounts 时发生错误: {e}"
        }

def collect_results_tool(srr_ids: List[str], **kwargs) -> dict:
    """
    收集分析结果
    """
    print(f"Tool 'collect_results_tool' called for SRRs {srr_ids}.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    
    # 收集各种结果文件
    results = {
        "fastp_results": [],
        "star_results": [],
        "featurecounts_results": []
    }
    
    for srr_id in srr_ids:
        # fastp 结果
        fastp_dir = data_dir / 'results' / 'fastp' / srr_id
        if fastp_dir.is_dir():
            results["fastp_results"].append({
                "srr_id": srr_id,
                "html_report": str(fastp_dir / f"{srr_id}.html"),
                "json_report": str(fastp_dir / f"{srr_id}.json")
            })
        
        # STAR 结果
        bam_dir = data_dir / 'results' / 'bam' / srr_id
        if bam_dir.is_dir():
            results["star_results"].append({
                "srr_id": srr_id,
                "bam_file": str(bam_dir / f"{srr_id}.bam"),
                "log_file": str(bam_dir / f"{srr_id}.log")
            })
    
    # featureCounts 结果
    featurecounts_dir = data_dir / 'results' / 'featurecounts'
    if featurecounts_dir.is_dir():
        results["featurecounts_results"] = {
            "counts_file": str(featurecounts_dir / "counts.txt"),
            "summary_file": str(featurecounts_dir / "counts.txt.summary")
        }
    
    return {
        "status": "success",
        "message": f"结果收集完成",
        "srr_ids": srr_ids,
        "results": results
    }

def generate_report_tool(srr_ids: List[str], **kwargs) -> dict:
    """
    生成分析报告
    """
    print(f"Tool 'generate_report_tool' called for SRR IDs: {srr_ids}")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    results_dir = data_dir / 'results'
    
    try:
        # 创建报告目录
        report_dir = results_dir / 'reports'
        report_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成HTML报告
        report_content = f"""
        <html>
        <head>
            <title>RNA-seq 分析报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 10px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .success {{ color: green; }}
                .error {{ color: red; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>RNA-seq 分析报告</h1>
                <p>分析时间: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>分析样本</h2>
                <ul>
        """
        
        for srr_id in srr_ids:
            report_content += f"<li>{srr_id}</li>\n"
        
        report_content += """
                </ul>
            </div>
            
            <div class="section">
                <h2>分析步骤</h2>
                <ul>
                    <li class="success">✓ 数据下载和验证</li>
                    <li class="success">✓ 质量控制 (fastp)</li>
                    <li class="success">✓ 序列比对 (STAR)</li>
                    <li class="success">✓ 基因定量 (featureCounts)</li>
                </ul>
            </div>
            
            <div class="section">
                <h2>结果文件</h2>
                <p>所有结果文件已保存在 <code>data/results/</code> 目录中。</p>
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
            "message": "分析报告生成完成",
            "report_file": str(report_file),
            "srr_ids": srr_ids
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

def react_plan_tool(analysis_type: str, **kwargs) -> dict:
    """
    React模式分析计划工具
    """
    print(f"Tool 'react_plan_tool' called for analysis type: {analysis_type}")
    
    if analysis_type == "rna_seq":
        plan = {
            "analysis_type": "rna_seq",
            "steps": [
                {"step": "check_environment", "tool": "check_environment_tool", "description": "检查分析环境"},
                {"step": "search_genome", "tool": "search_genome_tool", "description": "搜索参考基因组"},
                {"step": "download_genome", "tool": "download_genome_tool", "description": "下载基因组文件（如果需要）"},
                {"step": "build_index", "tool": "build_star_index_tool", "description": "构建STAR索引"},
                {"step": "search_fastq", "tool": "search_fastq_tool", "description": "搜索FASTQ文件"},
                {"step": "download_fastq", "tool": "download_fastq_tool", "description": "下载FASTQ文件（如果需要）"},
                {"step": "validate_fastq", "tool": "validate_fastq_tool", "description": "验证FASTQ文件"},
                {"step": "quality_control", "tool": "run_fastp_tool", "description": "质量控制"},
                {"step": "alignment", "tool": "run_star_align_tool", "description": "序列比对"},
                {"step": "quantification", "tool": "run_featurecounts_tool", "description": "基因定量"},
                {"step": "collect_results", "tool": "collect_results_tool", "description": "收集结果"},
                {"step": "generate_report", "tool": "generate_report_tool", "description": "生成报告"}
            ],
            "estimated_time": "2-4小时",
            "dependencies": ["conda", "wget", "fasterq-dump"]
        }
    else:
        plan = {
            "analysis_type": analysis_type,
            "error": f"不支持的分析类型: {analysis_type}"
        }
    
    return {
        "status": "success",
        "plan": plan,
        "react_mode": True
    }

def react_evaluate_tool(previous_result: dict, next_step: str, **kwargs) -> dict:
    """
    React模式结果评估工具
    """
    print(f"Tool 'react_evaluate_tool' called for step: {next_step}")
    
    # 分析前一步的结果
    status = previous_result.get("status", "unknown")
    message = previous_result.get("message", "")
    
    evaluation = {
        "previous_step_status": status,
        "previous_step_message": message,
        "next_step": next_step,
        "should_proceed": status == "success",
        "recommendation": ""
    }
    
    if status == "success":
        evaluation["recommendation"] = f"前一步执行成功，可以继续执行 {next_step}"
    elif status == "error":
        evaluation["recommendation"] = f"前一步执行失败: {message}，需要重新执行或跳过"
    elif status == "partial_success":
        evaluation["recommendation"] = f"前一步部分成功: {message}，可以继续但需要注意"
    else:
        evaluation["recommendation"] = f"前一步状态未知: {status}，建议重新检查"
    
    return {
        "status": "success",
        "evaluation": evaluation,
        "react_mode": True
    }

def react_summary_tool(completed_steps: List[str], failed_steps: List[str], **kwargs) -> dict:
    """
    React模式执行总结工具
    """
    print(f"Tool 'react_summary_tool' called.")
    
    summary = {
        "total_steps": len(completed_steps) + len(failed_steps),
        "completed_steps": completed_steps,
        "failed_steps": failed_steps,
        "success_rate": len(completed_steps) / (len(completed_steps) + len(failed_steps)) if (len(completed_steps) + len(failed_steps)) > 0 else 0,
        "overall_status": "success" if len(failed_steps) == 0 else "partial_success" if len(completed_steps) > 0 else "failed"
    }
    
    return {
        "status": "success",
        "summary": summary,
        "react_mode": True
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

def execute_rnaseq_workflow(srr_ids: List[str], genome_name: str, **kwargs) -> dict:
    """
    执行完整的RNA-seq分析工作流 - 通过Nextflow
    """
    print(f"Tool 'execute_rnaseq_workflow' called for SRRs {srr_ids} with genome '{genome_name}'.")
    
    # 构建完整的Nextflow命令
    project_root = _get_project_root()
    work_dir = project_root / 'data'
    genomes_config = _get_genomes_config()
    
    # 检查基因组配置
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        return {
            "status": "error",
            "message": f"基因组 '{genome_name}' 在配置中未找到"
        }
    
    try:
        # 检查输入文件
        fastq_dir = work_dir / 'fastq'
        fasta_path = project_root / genome_info['fasta']
        gtf_path = project_root / genome_info['gtf']
        
        # 构建Nextflow命令
        command = [
            "nextflow", "run", str(project_root / "nextflow" / "main.nf"),
            "--genome_name", genome_name,
            "--species", genome_info['species'],
            "--genome_version", genome_info['version'],
            "--outdir", str(work_dir / 'results'),
            "-resume"
        ]
        
        # 根据文件存在情况决定执行哪些步骤
        if not fasta_path.is_file() or not gtf_path.is_file():
            command.extend(["--run_genome_download", "true"])
            print(f"基因组文件不存在，将下载基因组: {genome_name}")
        else:
            command.extend(["--fasta", str(fasta_path), "--gtf", str(gtf_path)])
            print(f"使用现有基因组文件: {genome_name}")
        
        # 检查FASTQ文件
        fastq_files_exist = any((fastq_dir / f"{srr_id}_1.fastq.gz").is_file() for srr_id in srr_ids)
        if not fastq_files_exist:
            command.extend(["--run_download_fastq", "true", "--srr_ids", ",".join(srr_ids)])
            print(f"FASTQ文件不存在，将下载SRR: {srr_ids}")
        else:
            reads_pattern = str(fastq_dir / f"{{{','.join(srr_ids)}}}_*.fastq.gz")
            command.extend(["--reads", reads_pattern])
            print(f"使用现有FASTQ文件: {srr_ids}")
        
        # 总是执行分析步骤
        command.extend(["--run_analysis", "true"])
        
        # 设置环境变量
        env = os.environ.copy()
        env['NXF_HOME'] = str(work_dir / '.nextflow')
        
        print(f"执行Nextflow命令: {' '.join(command)}")
        
        # 执行命令
        result = subprocess.run(command, env=env, capture_output=True, text=True)
        
        if result.returncode == 0:
            return {
                "status": "success",
                "message": f"RNA-seq分析工作流完成 (通过Nextflow)",
                "srr_ids": srr_ids,
                "genome_name": genome_name,
                "output_dir": str(work_dir / 'results'),
                "stdout": result.stdout,
                "stderr": result.stderr
            }
        else:
            return {
                "status": "error",
                "message": f"RNA-seq分析工作流失败: {result.stderr}",
                "stdout": result.stdout,
                "stderr": result.stderr
            }
    except Exception as e:
        return {
            "status": "error",
            "message": f"执行RNA-seq工作流时发生错误: {e}"
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