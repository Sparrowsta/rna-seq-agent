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

    Args:
        params (Dict[str, Any]): 运行流程所需的参数字典。
                                 需要 'srr_list_path', 'fasta_path', 'gtf_path'。

    Returns:
        subprocess.Popen: 已启动的子进程对象。
    """
    srr_list_path = params.get("srr_list_path")
    fasta_path = params.get("fasta_path")
    gtf_path = params.get("gtf_path")

    if not all([srr_list_path, fasta_path, gtf_path]):
        logging.error(f"缺少必要的流程参数。srr_list, fasta, gtf 都必须提供。")
        raise ValueError("srr_list_path, fasta_path, and gtf_path are required.")
    
    if not os.path.exists(srr_list_path):
        logging.error(f"提供的SRR列表文件路径不存在: {srr_list_path}")
        raise ValueError("srr_list_path must be a valid file path.")

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
        "--srr_list", srr_list_path,
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
    # 用于直接测试此模块功能的示例代码
    print("正在测试 run_nextflow_pipeline 函数...")
    
    # 使用项目中已有的 SRR_list.txt 文件进行测试
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    test_srr_file = os.path.join(project_root, "data", "SRR_list.txt")
    
    print(f"使用测试文件: {test_srr_file}")

    if not os.path.exists(test_srr_file):
        print(f"错误：测试文件不存在！请确保 {test_srr_file} 存在。")
    else:
        try:
            test_params = {"srr_list_path": test_srr_file}
            proc = run_nextflow_pipeline(test_params)
            
            # 等待一小段时间，看看是否有即时错误
            time.sleep(2)
            
            # 读取输出
            # 注意：因为用了 -bg，主进程会很快退出，communicate 可能不会等到太多东西
            stdout, stderr = proc.communicate(timeout=10) 
            
            print("\n--- Nextflow STDOUT ---")
            print(stdout)
            print("\n--- Nextflow STDERR ---")
            print(stderr)
            
            if proc.returncode == 0:
                print("\n测试成功：Nextflow 命令已成功提交到后台运行。")
            else:
                print(f"\n测试失败：Nextflow 命令返回错误码 {proc.returncode}。")

        except Exception as e:
            print(f"测试过程中发生错误: {e}")