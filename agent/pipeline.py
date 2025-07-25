# agent/pipeline.py
import subprocess
import os
import logging
import time
from typing import Dict, Any

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_nextflow_pipeline(params: Dict[str, Any]) -> subprocess.Popen:
    """
    接收参数并异步启动一个 Nextflow 流程。
    此函数现在适配新的“责任分离”架构。

    Args:
        params (Dict[str, Any]): 运行流程所需的参数字典。
                                 需要 'reads_glob', 'fasta_path', 'gtf_path'。

    Returns:
        subprocess.Popen: 已启动的子进程对象。
    """
    reads_glob = params.get("reads_glob")
    fasta_path = params.get("fasta_path")
    gtf_path = params.get("gtf_path")

    if not all([reads_glob, fasta_path, gtf_path]):
        logging.error(f"缺少必要的流程参数。reads_glob, fasta_path, gtf_path 都必须提供。")
        raise ValueError("reads_glob, fasta_path, and gtf_path are required.")

    # 获取项目根目录，以便正确地找到 nextflow/main.nf
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    nextflow_script_path = os.path.join(project_root, "nextflow", "main.nf")

    # 构建 Nextflow 命令
    command = [
        "nextflow",
        "run",
        nextflow_script_path,
        "-profile",
        "standard",
        "--reads", reads_glob,
        "--fasta", fasta_path,
        "--gtf", gtf_path,
        "-resume", # 总是启用 resume 功能
        "-bg" # 让 Nextflow 在后台运行
    ]

    logging.info(f"准备执行 Nextflow 命令: {' '.join(command)}")

    # 使用 Popen 异步执行命令
    # stdout=subprocess.PIPE 和 stderr=subprocess.PIPE 可以捕获输出，以便后续查看
    process = subprocess.Popen(
        command,
        cwd=project_root,  # 在项目根目录执行
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True  # 让输出以文本形式返回
    )

    logging.info(f"Nextflow 流程已启动，进程 PID: {process.pid}")

    return process

if __name__ == '__main__':
    # 此模块的直接测试功能已过时，因为其依赖于动态的任务数据库和锁。
    # 请使用根目录下的 test_final_architecture.py 进行完整的端到端测试。
    pass