# agent/tools.py

import os
import time
import json
import threading
import subprocess
from pathlib import Path
from typing import List

# 为 genomes.json 创建一个专用的文件锁，防止并发写入冲突
genomes_json_lock = threading.Lock()


# 导入最底层的 pipeline 执行器
from agent.pipeline import run_nextflow_pipeline

def list_available_genomes() -> dict:
    """
    读取并返回 config/genomes.json 的内容，列出所有可用的基因组。
    """
    print("工具 'list_available_genomes' 被调用。")
    try:
        # 构造相对于项目根目录的绝对路径
        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'genomes.json')
        with open(config_path, 'r') as f:
            genomes = json.load(f)
        # 直接返回字典，server 端会处理序列化
        return genomes
    except FileNotFoundError:
        return {"error": "配置文件 'config/genomes.json' 未找到。"}
    except json.JSONDecodeError:
        return {"error": "配置文件 'config/genomes.json' 格式不正确。"}
    except Exception as e:
        return {"error": f"读取基因组配置时发生未知错误: {e}"}

def add_genome_to_config(genome_name: str, species: str, version: str, fasta_url: str, gtf_url: str) -> dict:
    """
    将一个新的基因组条目线程安全地添加到 'config/genomes.json'。
    这个函数只更新配置文件，不执行下载。
    """
    print(f"工具 'add_genome_to_config' 被调用，参数: {genome_name}")
    with genomes_json_lock:
        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'genomes.json')
        try:
            with open(config_path, 'r') as f:
                genomes_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            genomes_data = {}
        
        # 根据URL推断文件名，并构建相对路径
        fasta_filename = os.path.basename(fasta_url)
        gtf_filename = os.path.basename(gtf_url)
        relative_fasta_path = os.path.join('data', 'genomes', species, version, fasta_filename)
        relative_gtf_path = os.path.join('data', 'genomes', species, version, gtf_filename)

        genomes_data[genome_name] = {
            "species": species,
            "version": version,
            "fasta": relative_fasta_path,
            "gtf": relative_gtf_path,
            "fasta_url": fasta_url,
            "gtf_url": gtf_url
        }
        
        with open(config_path, 'w') as f:
            json.dump(genomes_data, f, indent=2)
        
        return {"status": "success", "message": f"基因组 '{genome_name}' 已成功添加到配置文件中。"}

def _perform_download(task_id: str, task_database: dict, db_lock: threading.Lock, genome_name: str, species: str, version: str, fasta_url: str, gtf_url: str):
    """在后台线程中执行实际的下载操作，并在完成后更新任务状态。"""
    print(f"后台任务 {task_id} 开始：为 '{genome_name}' 下载基因组文件...")
    base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'genomes', species, version)
    os.makedirs(base_path, exist_ok=True)

    fasta_local_path = os.path.join(base_path, os.path.basename(fasta_url))
    gtf_local_path = os.path.join(base_path, os.path.basename(gtf_url))

    try:
        print(f"正在下载 FASTA 文件从 {fasta_url} 到 {fasta_local_path}...")
        subprocess.run(['wget', '-O', fasta_local_path, fasta_url], check=True)
        print(f"FASTA 文件 '{genome_name}' 下载完成。")

        print(f"正在下载 GTF 文件从 {gtf_url} 到 {gtf_local_path}...")
        subprocess.run(['wget', '-O', gtf_local_path, gtf_url], check=True)
        print(f"GTF 文件 '{genome_name}' 下载完成。")

        # 下载成功，更新任务状态
        with db_lock:
            task_database[task_id]['status'] = 'completed'
            task_database[task_id]['end_time'] = time.time()

    except subprocess.CalledProcessError as e:
        print(f"错误：下载 '{genome_name}' 时 wget 命令失败: {e}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = str(e)
    except Exception as e:
        print(f"错误：执行下载任务 '{genome_name}' 时发生未知错误: {e}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = str(e)

def download_genome_files(genome_name: str, task_database: dict, db_lock: threading.Lock, task_id_counter: int) -> dict:
    """
    为一个在配置文件中已存在的基因组启动后台下载，并返回一个可追踪的任务ID。
    """
    print(f"工具 'download_genome_files' 被调用，准备下载 '{genome_name}'。")
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'genomes.json')
    
    try:
        with open(config_path, 'r') as f:
            genomes_data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {"status": "error", "message": "无法读取或解析基因组配置文件。"}

    genome_info = genomes_data.get(genome_name)
    if not genome_info:
        return {"status": "error", "message": f"在配置文件中找不到名为 '{genome_name}' 的基因组。"}

    # --- 新增：下载前检查文件是否存在 ---
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fasta_path = os.path.join(project_root, genome_info.get('fasta', ''))
    gtf_path = os.path.join(project_root, genome_info.get('gtf', ''))

    if os.path.exists(fasta_path) and os.path.exists(gtf_path):
        print(f"文件已存在，跳过下载 '{genome_name}'。")
        return {
            "status": "skipped",
            "message": f"基因组 '{genome_name}' 的文件已存在，跳过下载。"
        }
    # --- 检查结束 ---

    with db_lock:
        # 注意：这里的 task_id_counter 是一个 int，不能直接拼接，需要先转换
        task_id = f"download_{task_id_counter}"
        task_database[task_id] = {
            "type": "download",
            "status": "downloading",
            "start_time": time.time(),
            "details": {
                "genome_name": genome_name,
                "species": genome_info['species'],
                "version": genome_info['version'],
                "fasta_url": genome_info['fasta_url'],
                "gtf_url": genome_info['gtf_url']
            }
        }

    # 创建并启动一个后台线程来执行下载
    download_thread = threading.Thread(
        target=_perform_download,
        args=(
            task_id,
            task_database,
            db_lock,
            genome_name,
            genome_info['species'],
            genome_info['version'],
            genome_info['fasta_url'],
            genome_info['gtf_url']
        )
    )
    download_thread.daemon = True
    download_thread.start()
    
    # 注意：这里需要返回一个包含 updated_task_id_counter 的字典，以便 server 更新计数器
    return {
        "status": "success",
        "message": f"基因组 '{genome_name}' 的下载任务已在后台启动。",
        "task_id": task_id,
        "updated_task_id_counter": task_id_counter + 1
    }

def get_task_status(task_id: str, task_database: dict, db_lock: threading.Lock) -> dict:
    """
    查询指定任务ID的状态。
    对于正在运行的 pipeline 任务，此函数会检查底层进程的真实状态并更新数据库。
    """
    print(f"工具 'get_task_status' 被调用，查询任务: {task_id}")
    with db_lock:
        # 使用深拷贝复制一份，避免在检查过程中直接修改返回给用户的数据
        task_info = json.loads(json.dumps(task_database.get(task_id)))

    if not task_info:
        return {"status": "error", "message": f"找不到任务 '{task_id}'。"}

    # 对于 pipeline 任务，需要检查进程状态
    task_type = task_info.get("type")
    if task_type == "pipeline" and task_info.get("status") == "running":
        pid = task_info.get("details", {}).get("pid")
        if pid:
            try:
                # 非阻塞地检查进程状态
                pid_check, exit_status = os.waitpid(pid, os.WNOHANG)
                if pid_check != 0:
                    # 进程已结束，更新数据库
                    with db_lock:
                        # 再次检查以避免竞争条件
                        if task_id in task_database and task_database[task_id]['status'] == 'running':
                            exit_code = os.WEXITSTATUS(exit_status)
                            task_database[task_id]["status"] = "completed" if exit_code == 0 else "failed"
                            task_database[task_id]["details"]["exit_code"] = exit_code
                            task_database[task_id]["end_time"] = time.time()
                            print(f"任务 {task_id} (PID: {pid}) 已结束，退出码: {exit_code}。状态已更新。")
            except OSError:
                # 进程不存在，也视为任务失败
                with db_lock:
                    if task_id in task_database and task_database[task_id]['status'] == 'running':
                        task_database[task_id]["status"] = "failed"
                        task_database[task_id]["details"]["error"] = "Process not found."
                        task_database[task_id]["end_time"] = time.time()
                        print(f"任务 {task_id} (PID: {pid}) 的进程未找到。状态已更新为 failed。")

    # 在检查和潜在的更新后，返回最新的任务信息
    with db_lock:
        return task_database.get(task_id, {})

def _perform_rna_seq_pipeline(task_id: str, task_database: dict, db_lock: threading.Lock, genome_name: str, srr_list: str):
    """
    在后台线程中执行 RNA-seq 流程的完整操作，包括下载、压缩和运行 Nextflow。
    """
    # 注意：这个函数在自己的线程中运行，所以可以安全地执行阻塞操作。
    # 它通过闭包或参数捕获了所有需要的变量。

    # --- 1. 更新任务状态为 'running' 并解析参数 ---
    project_root = Path(__file__).parent.parent.resolve()
    config_path = project_root / 'config' / 'genomes.json'
    
    with db_lock:
        task_database[task_id]['status'] = 'running'
        task_database[task_id]['details']['message'] = '正在解析参数和基因组信息...'

    try:
        with open(config_path, 'r') as f:
            genomes_data = json.load(f)
        genome_info = genomes_data.get(genome_name)
        if not genome_info:
            raise ValueError(f"在配置文件中找不到名为 '{genome_name}' 的基因组。")
        
        fasta_path = str(project_root / genome_info['fasta'])
        gtf_path = str(project_root / genome_info['gtf'])
        srr_ids = srr_list.replace(",", " ").split()
        if not srr_ids:
            raise ValueError("SRR 列表不能为空。")

    except Exception as e:
        print(f"任务 {task_id} 失败: 参数解析错误。错误: {e}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = f"参数解析错误: {e}"
        return

    # --- 2. 检查、下载并压缩 FASTQ 文件 ---
    fastq_dir = project_root / 'data' / 'fastq'
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        with db_lock:
            task_database[task_id]['details']['message'] = '正在准备 FASTQ 文件...'
        
        for srr_id in srr_ids:
            with db_lock:
                 task_database[task_id]['details']['current_srr'] = srr_id
            
            existing_gz_files = list(fastq_dir.glob(f"{srr_id}*.fastq.gz"))
            if existing_gz_files:
                print(f"任务 {task_id}: 找到 '{srr_id}' 的本地文件，跳过下载。")
                continue

            print(f"任务 {task_id}: 未找到 '{srr_id}' 的本地文件，开始从 SRA 下载...")
            with db_lock:
                task_database[task_id]['details']['message'] = f"正在从 SRA 下载 {srr_id}..."

            dump_command = ['conda', 'run', '--no-capture-output', '-n', 'sra_env', 'fasterq-dump', '--split-files', '--outdir', str(fastq_dir), srr_id]
            dump_result = subprocess.run(dump_command, capture_output=True, text=True)
            if dump_result.returncode != 0:
                raise subprocess.CalledProcessError(dump_result.returncode, dump_command, dump_result.stdout, dump_result.stderr)
            
            print(f"任务 {task_id}: 成功下载 '{srr_id}'。")
            with db_lock:
                task_database[task_id]['details']['message'] = f"正在压缩 {srr_id}..."

            newly_dumped_files = list(fastq_dir.glob(f"{srr_id}*.fastq"))
            if not newly_dumped_files:
                print(f"警告 (任务 {task_id}): 为 '{srr_id}' 运行 fasterq-dump 后未找到输出的 .fastq 文件。")
                continue

            for fastq_file in newly_dumped_files:
                gzip_command = ['gzip', '-f', str(fastq_file)]
                gzip_result = subprocess.run(gzip_command, capture_output=True, text=True)
                if gzip_result.returncode != 0:
                    raise subprocess.CalledProcessError(gzip_result.returncode, gzip_command, gzip_result.stdout, gzip_result.stderr)
                print(f"任务 {task_id}: 成功压缩 {fastq_file.name}.gz。")

    except subprocess.CalledProcessError as e:
        error_message = f"处理 SRR ID '{srr_id}' 时命令 '{' '.join(e.cmd)}' 失败。\n返回码: {e.returncode}\n标准输出: {e.stdout}\n标准错误: {e.stderr}"
        print(f"任务 {task_id} 失败: {error_message}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = error_message
        return
    except Exception as e:
        print(f"任务 {task_id} 失败: 准备 FASTQ 文件时发生未知错误: {e}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = f"准备 FASTQ 文件时发生未知错误: {e}"
        return

    # --- 3. 构造 Nextflow 参数并执行 ---
    srr_glob_part = "{" + ",".join(srr_ids) + "}"
    reads_glob_pattern = str(fastq_dir / f"{srr_glob_part}_*.fastq.gz")
    
    try:
        with db_lock:
            task_database[task_id]['details']['message'] = '正在启动 Nextflow 流程...'
            task_database[task_id]['details']['reads_glob_pattern'] = reads_glob_pattern
        
        pipeline_params = {
            "reads_glob": reads_glob_pattern,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path
        }
        process = run_nextflow_pipeline(pipeline_params)
        
        print(f"任务 {task_id}: Nextflow 流程已启动，PID: {process.pid}")
        with db_lock:
            task_database[task_id]['details']['pid'] = process.pid
            task_database[task_id]['details']['message'] = 'Nextflow 流程正在运行。'
        
        # 这里不需要等待，get_task_status 会负责检查进程状态

    except Exception as e:
        print(f"任务 {task_id} 失败: 执行 Nextflow 流程时发生未知错误: {e}")
        with db_lock:
            task_database[task_id]['status'] = 'failed'
            task_database[task_id]['details']['error'] = f"执行 Nextflow 流程时发生未知错误: {e}"
        return

def run_rna_seq_pipeline(genome_name: str, srr_list: str, task_database: dict, db_lock: threading.Lock, task_id: str) -> dict:
    """
    根据指定的基因组和 SRR 列表，在后台启动 RNA-seq Nextflow 流程。
    此函数会立即返回一个任务 ID，实际工作在后台线程中进行。
    """
    print(f"工具 'run_rna_seq_pipeline' 被调用，准备启动任务 {task_id}。基因组: {genome_name}, SRR列表: {srr_list}")

    # --- 1. 在数据库中立即创建任务记录 ---
    with db_lock:
        task_database[task_id] = {
            "type": "pipeline",
            "status": "starting", # 初始状态
            "start_time": time.time(),
            "details": {
                "genome_name": genome_name,
                "srr_list": srr_list.replace(",", " ").split(),
                "message": "任务正在初始化..."
            }
        }

    # --- 2. 创建并启动后台线程 ---
    pipeline_thread = threading.Thread(
        target=_perform_rna_seq_pipeline,
        args=(task_id, task_database, db_lock, genome_name, srr_list)
    )
    pipeline_thread.daemon = True
    pipeline_thread.start()

    # --- 3. 立即返回成功信息 ---
    return {
        "status": "success",
        "message": f"成功启动 RNA-seq 流程！任务 ID: {task_id}",
        "task_id": task_id
    }

def list_files(path: str = ".") -> dict:
    """
    列出指定路径下的文件和目录。
    为了安全，路径被限制在 'data' 目录内。
    例如，要查看 data/genomes/human 目录，请使用 path="genomes/human"。
    """
    print(f"工具 'list_files' 被调用，路径: {path}")
    
    # 安全限制：将根目录固定为项目下的 'data' 目录
    project_root = Path(__file__).parent.parent.resolve()
    data_root = project_root / 'data'
    
    # 构造并清理目标路径
    # 使用 resolve() 来处理 ".." 等路径，并确保它是一个绝对路径
    target_path = (data_root / path).resolve()
    
    # 再次检查，确保最终路径在 data_root 内，防止路径逃逸
    if not target_path.is_relative_to(data_root.resolve()):
         return {"error": "访问被拒绝。只能查看 'data' 目录下的内容。"}
        
    try:
        if not target_path.is_dir():
            return {"error": f"路径 '{path}' 不是一个有效的目录。"}
            
        items = [item.name for item in target_path.iterdir()]
        return {"path": str(target_path.relative_to(data_root)), "files": items}
    except FileNotFoundError:
        return {"error": f"路径 '{path}' 未找到。"}
    except Exception as e:
        return {"error": f"列出文件时发生错误: {e}"}

def unsupported_request(user_request: str) -> dict:
    """
    当用户的请求与任何其他可用工具的功能都不匹配时，必须调用此工具。
    它会向用户返回一个标准化的消息，说明请求无法处理。
    """
    print(f"工具 'unsupported_request' 被调用，原始请求: {user_request}")
    # 在未来，这里可以加入更复杂的逻辑，比如记录无法处理的请求用于分析
    return {
        "status": "error",
        "type": "unsupported_request",
        "message": "抱歉，我无法处理您的请求。我的能力目前仅限于运行生物信息学流程和管理相关数据。请尝试提出与以下功能相关的请求：运行分析、查询任务状态、列出可用基因组、添加或下载基因组、查看文件列表。"
    }

# 你可以在这里直接测试这个工具
if __name__ == '__main__':
    # 这个 __main__ 块需要更新以反映新的函数签名和逻辑
    # 由于它依赖于任务数据库等，直接运行变得复杂，因此暂时将其留空。
    pass