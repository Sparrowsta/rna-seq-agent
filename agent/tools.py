# agent/tools.py - v5.2 Architecture
import os
import time
import json
import threading
from pathlib import Path
from typing import Dict, Any

# --- New Imports ---
from agent.executor import execute_with_logging
from agent.log_parser import parse_log_line

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

        genome_info = plan['genome_info']
        command.extend(["--fasta", str(project_root / genome_info['fasta'])])
        command.extend(["--gtf", str(project_root / genome_info['gtf'])])
        command.extend(["--species", genome_info['species']])
        command.extend(["--genome_version", genome_info['version']])

        if 'download_fastq' in plan['steps_to_execute']:
            command.extend(["--srr_ids", ",".join(plan['srr_ids'])])
        else:
            srr_glob = "{" + ",".join(plan['srr_ids']) + "}"
            reads_pattern = str(work_dir / 'fastq' / f"{srr_glob}_*_{'{1,2}'}.fastq.gz")
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
        process = execute_with_logging(command, log_file, work_dir)
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
        "is_executable": True,
        "genome_info": None
    }

    # 1. Check genome dependencies
    genomes_config = _get_genomes_config()
    genome_info = genomes_config.get(genome_name)
    if not genome_info:
        plan['messages'].append(f"❌ **错误**: 在配置中找不到基因组 '{genome_name}'。")
        plan['is_executable'] = False
        return {"status": "success", "plan": plan}
    
    plan['genome_info'] = genome_info
    index_path = project_root / genome_info['fasta']
    index_path = index_path.parent / 'star_index'
    
    if not index_path.is_dir():
        plan['steps_to_execute'].append('genome_prep')
        plan['messages'].append(f"ℹ️ **计划**: 将为 '{genome_name}' 新建 STAR 索引。")
    else:
        plan['messages'].append(f"✅ **检查通过**: '{genome_name}' 的 STAR 索引已存在。")

    # 2. Check FASTQ file dependencies
    missing_srr_files = [srr for srr in plan['srr_ids'] if not (data_dir / 'fastq' / f"{srr}_1.fastq.gz").exists()]
    
    if missing_srr_files:
        plan['steps_to_execute'].append('download_fastq')
        plan['messages'].append(f"ℹ️ **计划**: 将下载 {len(missing_srr_files)} 个缺失的样本文件: {', '.join(missing_srr_files)}。")
    else:
        plan['messages'].append("✅ **检查通过**: 所有需要的 FASTQ 文件均已本地存在。")

    # 3. Analysis is always the final step
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
        return json.loads(json.dumps(task_info))

def list_available_genomes() -> dict:
    """Reads and returns the content of config/genomes.json."""
    print("Tool 'list_available_genomes' called.")
    genomes_data = _get_genomes_config()
    if not genomes_data:
        return {"error": "Could not read or find 'config/genomes.json'."}
    return {"status": "success", "genomes": genomes_data}