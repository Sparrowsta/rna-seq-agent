# agent/server.py

import json
import time
import asyncio
import os
import tempfile
import re
from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import List, AsyncGenerator, Dict, Any
from threading import Lock
from fastapi.middleware.cors import CORSMiddleware

# 导入我们新创建的 pipeline 模块和 tools 模块
from agent.pipeline import run_nextflow_pipeline
from agent.tools import run_rna_seq_pipeline, get_task_status

# --- 1. 定义数据模型 ---
class StandardMessage(BaseModel):
    role: str
    content: str

class ChatInput(BaseModel):
    messages: List[StandardMessage]

class PipelineRunInput(BaseModel):
    srr_list: str

# --- 2. 创建 FastAPI 应用实例 ---
app = FastAPI(
    title="从零构建的流式 Agent 服务器",
    description="这是一个支持流式响应、符合OpenAI接口标准的Agent服务器。",
    version="0.3.0", # 版本更新：增加了任务跟踪
)

# --- 3. 添加 CORS 中间件 ---
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- 3.5. 全局任务数据库 (内存中) ---
# 使用字典来存储任务状态。在生产环境中，这应该被替换为真正的数据库。
TASK_DATABASE: Dict[str, Dict[str, Any]] = {}
# 使用锁来确保对 TASK_DATABASE 的线程安全访问
db_lock = Lock()
# 任务 ID 计数器
task_id_counter = 1


# --- 4. 定义 API 端点 ---

@app.get("/")
def read_root():
    """根路径，用于检查服务器是否在线。"""
    return {"message": "你好，流式 Agent 服务器正在运行！"}

@app.get("/models") 
async def list_models():
    """
    响应客户端获取模型列表的请求。
    我们使用自己定义的、清晰的 model_id。
    """
    model_id = "rna-seq-agent-v1"
    return {
        "object": "list",
        "data": [
            {
                "id": model_id,
                "object": "model",
                "created": int(time.time()),
                "owned_by": "user",
            }
        ],
    }

async def stream_agent_response(chat_input: ChatInput) -> AsyncGenerator[str, None]:
    """
    一个真正的 Agent 响应生成器。
    它会分析用户输入，并决定是调用工具还是返回普通信息。
    """
    model_id = "rna-seq-agent-v1"
    user_message = chat_input.messages[-1].content if chat_input.messages else ""

    # --- Agent 核心逻辑：路由和工具调用 ---
    
    # 规则1：如果用户想要运行流程
    # 我们使用非常简单的正则匹配来查找 SRR ID 和关键词
    srr_pattern = re.compile(r"(SRR\d{7,})", re.IGNORECASE)
    run_keywords = ["run", "运行", "执行", "处理"]
    
    srr_matches = srr_pattern.findall(user_message)
    contains_run_keyword = any(keyword in user_message for keyword in run_keywords)

    # 规则2：如果用户想要查询任务状态
    task_pattern = re.compile(r"(task_\d+)", re.IGNORECASE)
    query_keywords = ["查询", "状态", "status", "检查", "check"]
    
    task_matches = task_pattern.findall(user_message)
    contains_query_keyword = any(keyword in user_message for keyword in query_keywords)

    response_text = ""
    if srr_matches and contains_run_keyword:
        # 如果找到了 SRR ID 和运行关键词，就调用工具
        srr_list_str = ", ".join(srr_matches)
        print(f"检测到运行意图，提取到 SRR IDs: {srr_list_str}")
        
        tool_result = run_rna_seq_pipeline(srr_list_str)
        
        if tool_result.get("status") == "success":
            global task_id_counter
            with db_lock:
                task_id = f"task_{task_id_counter}"
                task_id_counter += 1
                
                TASK_DATABASE[task_id] = {
                    "pid": tool_result.get("pid"),
                    "srr_list_file": tool_result.get("srr_list_file"),
                    "status": "running", # 初始状态
                    "start_time": time.time(),
                    "srr_list": srr_list_str
                }
            
            response_text = f"任务 '{task_id}' 已成功启动。进程 PID: {tool_result.get('pid')}。"
        else:
            # 如果工具返回错误
            error_message = tool_result.get("message", "未知工具错误")
            response_text = f"启动流程失败: {error_message}"

    elif task_matches and contains_query_keyword:
        # 如果找到了 task_id 和查询关键词，就查询状态
        task_id_to_query = task_matches[0]
        print(f"检测到查询意图，查询任务 ID: {task_id_to_query}")
        
        with db_lock:
            task_info = TASK_DATABASE.get(task_id_to_query)

        if not task_info:
            response_text = f"错误：找不到任务 '{task_id_to_query}'。"
        else:
            pid = task_info.get("pid")
            current_status = task_info.get("status", "unknown")
            
            # 使用 os.waitpid 进行更精确、更健壮的状态检查
            process_status_text = "未知"
            if current_status == "running" and pid:
                try:
                    pid_check, exit_status = os.waitpid(pid, os.WNOHANG)
                    if pid_check == 0:
                        process_status_text = "仍在运行"
                    else:
                        exit_code = os.WEXITSTATUS(exit_status)
                        process_status_text = f"已结束 (退出码: {exit_code})"
                        with db_lock:
                            TASK_DATABASE[task_id_to_query]["status"] = "finished"
                            TASK_DATABASE[task_id_to_query]["exit_code"] = exit_code
                except OSError:
                    process_status_text = "已结束 (进程不存在)"
                    with db_lock:
                        if TASK_DATABASE.get(task_id_to_query):
                            TASK_DATABASE[task_id_to_query]["status"] = "finished"

            elif current_status == "finished":
                exit_code = task_info.get('exit_code', 'N/A')
                process_status_text = f"已结束 (退出码: {exit_code})"
            
            response_text = f"任务 '{task_id_to_query}' 的状态为: {process_status_text} (PID: {pid})。"

    else:
        # 否则，返回默认的帮助信息
        response_text = "你好！我是 RNA-seq 流程助手。你可以通过发送'运行 SRR1234567'或'查询 task_1'来与我交互。"

    # --- 将 Agent 的响应文本转换为流式输出 ---
    
    # 将完整响应文本逐字发送
    for char in response_text:
        response_chunk = {
            "id": f"chatcmpl-{int(time.time())}",
            "object": "chat.completion.chunk",
            "created": int(time.time()),
            "model": model_id,
            "choices": [{
                "index": 0,
                "delta": {"content": char},
                "finish_reason": None
            }]
        }
        sse_formatted_chunk = f"data: {json.dumps(response_chunk)}\n\n"
        yield sse_formatted_chunk
        await asyncio.sleep(0.02) # 模拟打字效果

    # 发送结束标志
    final_chunk = {
        "id": f"chatcmpl-{int(time.time())}",
        "object": "chat.completion.chunk",
        "created": int(time.time()),
        "model": model_id,
        "choices": [{
            "index": 0,
            "delta": {},
            "finish_reason": "stop"
        }]
    }
    yield f"data: {json.dumps(final_chunk)}\n\n"
    yield "data: [DONE]\n\n"

@app.post("/chat/completions")
async def chat_completions(chat_input: ChatInput):
    """
    核心聊天接口，现在它会调用真正的 Agent 逻辑。
    """
    return StreamingResponse(stream_agent_response(chat_input), media_type="text/event-stream")

@app.get("/tasks/{task_id}")
async def get_task_status(task_id: str):
    """
    根据任务 ID 查询任务的状态和详细信息。
    """
    with db_lock:
        task = TASK_DATABASE.get(task_id)

    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    # 使用 os.waitpid 进行更精确、更健壮的状态检查
    pid = task.get("pid")
    if pid and task.get("status") == "running":
        try:
            # os.WNOHANG 使其成为非阻塞调用
            pid_check, exit_status = os.waitpid(pid, os.WNOHANG)
            if pid_check == 0:
                # 子进程仍在运行
                task["process_status"] = "running"
            else:
                # 子进程已结束
                task["process_status"] = "finished"
                task["exit_code"] = os.WEXITSTATUS(exit_status) # 获取退出码
                # 更新数据库中的状态
                with db_lock:
                    TASK_DATABASE[task_id]["status"] = "finished"
        except OSError:
            # 如果进程不存在（例如，已经被其他方式清理），则也视为已结束
            task["process_status"] = "finished_or_not_found"
            with db_lock:
                if TASK_DATABASE.get(task_id):
                    TASK_DATABASE[task_id]["status"] = "finished"
    elif task.get("status") == "finished":
        task["process_status"] = "finished"
    
    return task

# @app.post("/run_pipeline")
# async def run_pipeline_endpoint(run_input: PipelineRunInput):
#     """
#     [已弃用] 用于启动 Nextflow 流程的 API 端点。
#     这个端点的业务逻辑已经全部移动到 agent/tools.py 的 run_rna_seq_pipeline 函数中，
#     以解决死锁问题并优化代码结构。Agent 现在直接调用工具函数，不再需要这个 API。
#     保留此代码作为参考。
#     """
#     # 创建一个临时目录来存放 SRR 列表文件，如果尚不存在
#     temp_dir = os.path.join(os.path.dirname(__file__), "temp_srr_lists")
#     os.makedirs(temp_dir, exist_ok=True)
#
#     try:
#         # 将逗号或空格分隔的 SRR 字符串转换为列表
#         srr_ids = run_input.srr_list.replace(",", " ").split()
#         if not srr_ids:
#             raise HTTPException(status_code=400, detail="SRR list cannot be empty.")
#
#         # 创建一个带时间戳的、唯一的临时文件名
#         timestamp = int(time.time())
#         temp_file_path = os.path.join(temp_dir, f"srr_list_{timestamp}.txt")
#
#         # 将 SRR ID 写入临时文件，每行一个
#         with open(temp_file_path, "w") as f:
#             for srr_id in srr_ids:
#                 f.write(f"{srr_id}\n")
#
#         # 调用 pipeline 模块中的函数，传递临时文件的路径
#         process = run_nextflow_pipeline({"srr_list_path": temp_file_path})
#
#         return {
#             "message": "Nextflow pipeline started successfully.",
#             "pid": process.pid,
#             "srr_list_file": temp_file_path # 返回使用的临时文件路径，便于调试
#         }
#     except ValueError as e:
#         raise HTTPException(status_code=400, detail=str(e))
#     except Exception as e:
#         # 其他未知错误，返回 500 错误
#         raise HTTPException(status_code=500, detail=f"An unexpected error occurred: {e}")

# --- 5. 启动服务器 ---
if __name__ == "__main__":
    import uvicorn
    print("正在启动流式 Agent 服务器 v0.3.0 ，访问 http://localhost:8001")
    uvicorn.run(app, host="0.0.0.0", port=8001)