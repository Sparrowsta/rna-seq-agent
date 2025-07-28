# agent/tools.py - v5.2 Architecture (with all tools implemented)
import os
import time
import json
import threading
from pathlib import Path
from typing import Dict, Any, List

# --- New Imports ---
from agent.executor import execute_with_logging
from agent.log_parser import parse_log_line

# 移除进度条模块导入
# from .progress_bar import format_progress_update, get_progress_summary

# --- Module-level lock for config file operations ---
_config_lock = threading.Lock()

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
        
        for step in plan['steps_to_execute']:
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
        # This makes the system robust against LLM truncating the plan['steps_to_execute']
        fasta_path = project_root / genome_info['fasta']
        gtf_path = project_root / genome_info['gtf']
        files_exist = fasta_path.is_file() and gtf_path.is_file()
        
        # Only add genome_download if files don't exist AND it's not already in the plan
        if not files_exist and 'genome_download' not in plan['steps_to_execute']:
            print(f"[{task_id}] WARNING: Genome files don't exist but 'genome_download' not in plan. Adding it.")
            plan['steps_to_execute'].insert(0, 'genome_download')  # Add at beginning
        elif files_exist and 'genome_download' in plan['steps_to_execute']:
            print(f"[{task_id}] INFO: Genome files exist, removing 'genome_download' from plan.")
            plan['steps_to_execute'].remove('genome_download')  # Remove if files exist

        # Only provide FASTA/GTF paths if we are NOT downloading them,
        # meaning they must already exist.
        if 'genome_download' not in plan['steps_to_execute']:
            command.extend(["--fasta", str(project_root / genome_info['fasta'])])
            command.extend(["--gtf", str(project_root / genome_info['gtf'])])

        # Provide SRR IDs for download, or a glob pattern for existing reads.
        if 'download_fastq' in plan['steps_to_execute']:
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
    Performs a pre-flight check for a given set of SRR IDs and a genome.
    It checks for the existence of required files (STAR index, FASTQ files)
    and returns a detailed, human-readable execution plan.
    This tool does NOT execute any long-running tasks.
    """
    print(f"Tool 'plan_analysis_task' called for genome '{genome_name}' and SRRs '{srr_ids}'.")
    project_root = _get_project_root()
    data_dir = project_root / 'data'
    
    plan = {
        "srr_ids": sorted(list(set(srr_ids.replace(",", " ").split()))),
        "genome_name": genome_name,
        "steps_to_execute": [],
        "messages": [],
        "is_executable": False, # Default to not executable
        "genome_info": None
    }

    # --- Step 1: Validate Genome Configuration ---
    genomes_config = _get_genomes_config()
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        plan['messages'].append(f"❌ **错误**: 在配置中找不到基因组 '{genome_name}'。请从可用基因组列表中选择。")
        return {"status": "success", "plan": plan}
    # CRITICAL: Store the original, complete genome_info block.
    # This ensures that all keys (fasta, gtf, fasta_url, gtf_url) are always
    # present in the plan, preventing downstream KeyErrors.
    plan['genome_info'] = genome_info.copy()

    # --- Step 2: Check for Genome Source Files and Plan Download if Needed ---
    fasta_path = project_root / plan['genome_info']['fasta']
    gtf_path = project_root / plan['genome_info']['gtf']
    
    source_files_exist = fasta_path.is_file() and gtf_path.is_file()
    if not source_files_exist:
        plan['steps_to_execute'].append('genome_download')
        plan['messages'].append(f"ℹ️ **计划**: 将首先从 UCSC 下载 '{genome_name}' 的基因组源文件。")
    else:
        plan['messages'].append(f"✅ **基因组源文件**: '{genome_name}' 的 FASTA 和 GTF 文件已找到。")
        # Don't add genome_download step if files already exist

    # --- If we reach here, the plan is potentially executable ---
    plan['is_executable'] = True

    # --- Step 3: Check for STAR Index ---
    index_path = fasta_path.parent / 'star_index'
    if not (index_path / 'SA').is_file():
        plan['steps_to_execute'].append('genome_prep')
        plan['messages'].append(f"ℹ️ **计划**: 将为 '{genome_name}' 新建 STAR 索引。")
    else:
        plan['messages'].append(f"✅ **STAR 索引**: 已找到且完整。")

    # --- Step 4: Check for FASTQ Files ---
    missing_srr_files = [srr for srr in plan['srr_ids'] if not (data_dir / 'fastq' / f"{srr}_1.fastq.gz").exists()]
    if missing_srr_files:
        plan['steps_to_execute'].append('download_fastq')
        plan['messages'].append(f"ℹ️ **计划**: 将下载 {len(missing_srr_files)} 个缺失的样本文件: {', '.join(missing_srr_files)}。")
    else:
        plan['messages'].append("✅ **样本文件**: 所有需要的 FASTQ 文件均已本地存在。")

    # --- Step 5: Add Core Analysis Step ---
    plan['steps_to_execute'].append('analysis')
    plan['messages'].append("➡️ **最终步骤**: 将对所有样本执行核心分析 (QC, 比对, 定量)。")

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
        
        steps = []
        if 'genome_download' in plan['steps_to_execute']:
            steps.append({"name": "WGET_GENOME_FILES", "status": "pending", "logs": []})
        if 'genome_prep' in plan['steps_to_execute']:
            steps.append({"name": "STAR_GENOME_GENERATE", "status": "pending", "logs": []})
        if 'download_fastq' in plan['steps_to_execute']:
            steps.append({"name": "FASTERQ_DUMP", "status": "pending", "logs": []})
        if 'analysis' in plan['steps_to_execute']:
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
        "status": "success", "message": f"Task {task_id} has been successfully launched.",
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
        "steps_to_execute": ['genome_download', 'genome_prep'],
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
    When the user's request does not match any other available tools,
    this tool must be called. It returns a standardized message to the user
    indicating that the request cannot be handled.
    """
    print(f"Tool 'unsupported_request' called with original request: {user_request}")
    return {
        "status": "error",
        "type": "unsupported_request",
        "message": "I'm sorry, I cannot handle your request. My capabilities are currently focused on running bioinformatics pipelines and managing related data. Please try asking something related to: running an analysis, checking task status, or listing available genomes."
    }