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

def _perform_download(genome_name: str, species: str, version: str, fasta_url: str, gtf_url: str):
    """在后台线程中执行实际的下载操作。"""
    print(f"后台任务开始：为 '{genome_name}' 下载基因组文件...")
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
    except subprocess.CalledProcessError as e:
        print(f"错误：下载 '{genome_name}' 时 wget 命令失败: {e}")
    except Exception as e:
        print(f"错误：执行下载任务 '{genome_name}' 时发生未知错误: {e}")

def download_genome_files(genome_name: str) -> dict:
    """
    为一个在配置文件中已存在的基因组启动后台下载。
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

    # 创建并启动一个后台线程来执行下载
    download_thread = threading.Thread(
        target=_perform_download,
        args=(
            genome_name,
            genome_info['species'],
            genome_info['version'],
            genome_info['fasta_url'],
            genome_info['gtf_url']
        )
    )
    download_thread.daemon = True
    download_thread.start()
    
    return {
        "status": "success",
        "message": f"基因组 '{genome_name}' 的下载任务已在后台启动。"
    }

def run_rna_seq_pipeline(srr_list: str) -> dict:
    """
    调用 RNA-seq Nextflow 流程的工具 (厨师)。
    这个函数现在封装了所有的业务逻辑：创建临时文件并直接调用 pipeline 执行器。

    Args:
        srr_list: 一个包含 SRR 运行编号的字符串，可以由逗号或空格分隔。
                  例如: "SRR17469059, SRR17469060"

    Returns:
        一个包含操作结果的字典。
        成功: {'status': 'success', 'pid': ..., 'srr_list_file': '...'}
        失败: {'status': 'error', 'message': '...'}
    """
    print(f"工具 'run_rna_seq_pipeline' 被调用，参数: {srr_list}")
    
    # --- 核心业务逻辑 (准备食材) ---
    
    # 创建一个临时目录来存放 SRR 列表文件，如果尚不存在
    # 使用 __file__ 来确保路径总是相对于当前文件位置，更健壮
    temp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "temp_srr_lists")
    os.makedirs(temp_dir, exist_ok=True)

    try:
        # 将逗号或空格分隔的 SRR 字符串转换为列表
        srr_ids = srr_list.replace(",", " ").split()
        if not srr_ids:
            return {"status": "error", "message": "SRR 列表不能为空。"}

        # 创建一个带时间戳的、唯一的临时文件名
        timestamp = int(time.time())
        temp_file_path = os.path.join(temp_dir, f"srr_list_{timestamp}.txt")

        # 将 SRR ID 写入临时文件，每行一个
        with open(temp_file_path, "w") as f:
            for srr_id in srr_ids:
                f.write(f"{srr_id}\n")

        # --- 调用底层驱动 (使用炉灶) ---
        process = run_nextflow_pipeline({"srr_list_path": temp_file_path})
        
        return {
            "status": "success",
            "message": "成功启动 RNA-seq 流程！",
            "pid": process.pid,
            "srr_list_file": temp_file_path
        }

    except Exception as e:
        # 捕获所有可能的异常，例如文件系统权限问题
        return {"status": "error", "message": f"执行流程时发生未知错误: {e}"}

# 你可以在这里直接测试这个工具
if __name__ == '__main__':
    # 示例：如何调用这个工具
    # 这个测试现在不再需要 server 运行，可以直接执行
    test_srr_ids = "SRR17469059,SRR17469061"
    result_message = run_rna_seq_pipeline(test_srr_ids)
    print(result_message)