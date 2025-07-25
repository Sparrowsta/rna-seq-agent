# agent/tools.py

import os
import time
import json
import threading
import subprocess

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

def run_rna_seq_pipeline(genome_name: str, srr_list: str, task_database: dict, db_lock: threading.Lock, task_id: str) -> dict:
    """
    根据指定的基因组和 SRR 列表，调用 RNA-seq Nextflow 流程，并在任务数据库中创建一个条目。
    """
    print(f"工具 'run_rna_seq_pipeline' 被调用，基因组: {genome_name}, SRR列表: {srr_list}")

    # --- 1. 查找基因组文件路径 ---
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config', 'genomes.json')
    try:
        with open(config_path, 'r') as f:
            genomes_data = json.load(f)
        genome_info = genomes_data.get(genome_name)
        if not genome_info:
            return {"status": "error", "message": f"在配置文件中找不到名为 '{genome_name}' 的基因组。"}
        
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        fasta_path = os.path.join(project_root, genome_info['fasta'])
        gtf_path = os.path.join(project_root, genome_info['gtf'])

    except (FileNotFoundError, json.JSONDecodeError):
        return {"status": "error", "message": "无法读取或解析基因组配置文件。"}

    # --- 2. 准备 SRR 列表文件 ---
    temp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp_srr_lists")
    os.makedirs(temp_dir, exist_ok=True)
    try:
        srr_ids = srr_list.replace(",", " ").split()
        if not srr_ids:
            return {"status": "error", "message": "SRR 列表不能为空。"}

        timestamp = int(time.time())
        temp_file_path = os.path.join(temp_dir, f"srr_list_{timestamp}.txt")

        with open(temp_file_path, "w") as f:
            for srr_id in srr_ids:
                f.write(f"{srr_id}\n")

        # --- 3. 调用底层流程执行器 ---
        pipeline_params = {
            "srr_list_path": temp_file_path,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path
        }
        process = run_nextflow_pipeline(pipeline_params)
        
        # --- 4. 在数据库中创建任务记录 ---
        with db_lock:
            task_database[task_id] = {
                "type": "pipeline",
                "status": "running",
                "start_time": time.time(),
                "details": {
                    "pid": process.pid,
                    "genome_name": genome_name,
                    "srr_list": srr_ids,
                    "srr_list_file": temp_file_path
                }
            }
        
        return {
            "status": "success",
            "message": f"成功启动 RNA-seq 流程！任务 ID: {task_id}",
            "task_id": task_id
        }

    except Exception as e:
        return {"status": "error", "message": f"执行流程时发生未知错误: {e}"}

def list_files(path: str = ".") -> dict:
    """
    列出指定路径下的文件和目录。
    为了安全，路径被限制在 'data' 目录内。
    例如，要查看 data/genomes/human 目录，请使用 path="genomes/human"。
    """
    print(f"工具 'list_files' 被调用，路径: {path}")
    
    # 安全限制：将根目录固定为项目下的 'data' 目录
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_root = os.path.join(project_root, 'data')
    
    # 构造并清理目标路径


    target_path = os.path.join(data_root, path)
    
    # 再次检查，确保最终路径在 data_root 内，防止 ".." 等路径逃逸
    if not os.path.abspath(target_path).startswith(os.path.abspath(data_root)):
        return {"error": "访问被拒绝。只能查看 'data' 目录下的内容。"}
        
    try:
        if not os.path.isdir(target_path):
            return {"error": f"路径 '{path}' 不是一个有效的目录。"}
            
        items = os.listdir(target_path)
        return {"path": path, "files": items}
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
    # 示例：如何调用这个工具
    # 这个测试现在不再需要 server 运行，可以直接执行
    test_srr_ids = "SRR17469059,SRR17469061"
    result_message = run_rna_seq_pipeline(test_srr_ids)
    print(result_message)