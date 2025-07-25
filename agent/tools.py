# agent/tools.py

import os
import time

# 导入最底层的 pipeline 执行器
from agent.pipeline import run_nextflow_pipeline

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