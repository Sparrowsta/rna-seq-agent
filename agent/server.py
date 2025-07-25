# agent/server.py

import json
import time
import asyncio
import os
import tempfile
import re
import json
from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import List, AsyncGenerator, Dict, Any
from threading import Lock
from fastapi.middleware.cors import CORSMiddleware

# --- 新增: LLM 和环境相关导入 ---
import openai
from dotenv import load_dotenv

# --- 新增: 从我们的模块导入 ---
from agent.prompt import SYSTEM_PROMPT, TOOLS
import agent.tools as tool_module # 导入整个模块以便于函数查找

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
    title="LLM 驱动的流式 Agent 服务器",
    description="这是一个由大型语言模型驱动、支持工具调用、符合OpenAI接口标准的Agent服务器。",
    version="1.0.0", # 版本跃迁
)

# --- 3. 配置和初始化 ---
# 加载 .env 文件中的环境变量
# 这应该在访问任何环境变量之前完成
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), '..', 'config', '.env'))

# 初始化 OpenAI 客户端
# 从环境变量中获取配置
api_key = os.getenv("OPENAI_API_KEY")
base_url = os.getenv("OPENAI_API_BASE")
model_name = os.getenv("OPENAI_MODEL_NAME", "default-model") # 提供一个默认值

if not api_key or not base_url:
    raise ValueError("请在 config/.env 文件中设置 OPENAI_API_KEY 和 OPENAI_API_BASE")

client = openai.OpenAI(
    api_key=api_key,
    base_url=base_url,
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
    现在模型列表直接从环境变量中读取。
    """
    return {
        "object": "list",
        "data": [
            {
                "id": model_name,
                "object": "model",
                "created": int(time.time()),
                "owned_by": "system",
            }
        ],
    }

async def stream_agent_response(chat_input: ChatInput) -> AsyncGenerator[str, None]:
    """
    一个由 LLM 驱动的、支持工具调用的 Agent 响应生成器。
    """
    # 1. 准备发送给 LLM 的消息，确保我们的系统提示是唯一的
    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    # 过滤掉任何可能从上游传入的 system 消息，只保留 user 和 assistant 的消息
    for msg in chat_input.messages:
        if msg.role != "system":
            messages.append({"role": msg.role, "content": msg.content})

    # 2. 第一次调用 LLM，让它决定是否使用工具
    print(f"--- 第一次调用 LLM (模型: {model_name}) ---")
    print(f"发送的消息: {messages}")
    print(f"可用的工具: {TOOLS}")
    
    try:
        response = client.chat.completions.create(
            model=model_name,
            messages=messages,
            tools=TOOLS,
            tool_choice="required",
        )
    except Exception as e:
        print(f"调用 LLM API 时发生错误: {e}")
        yield f"data: {json.dumps({'error': str(e)})}\n\n"
        yield "data: [DONE]\n\n"
        return

    response_message = response.choices[0].message
    tool_calls = response_message.tool_calls

    # 3. 检查 LLM 是否决定调用工具
    if tool_calls:
        print(f"--- LLM 决定调用工具: {tool_calls} ---")
        messages.append(response_message)  # 将 assistant 的回复（包括工具调用请求）添加到历史记录中

        # 4. 执行所有工具调用
        for tool_call in tool_calls:
            function_name = tool_call.function.name
            
            try:
                function_args = json.loads(tool_call.function.arguments)
            except json.JSONDecodeError:
                function_response = f"错误: LLM 返回了无效的 JSON 参数: {tool_call.function.arguments}"
                messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                continue

            print(f"执行工具: {function_name}，参数: {function_args}")

            available_tools = {
                "run_rna_seq_pipeline": tool_module.run_rna_seq_pipeline,
                "list_available_genomes": tool_module.list_available_genomes,
                "add_genome_to_config": tool_module.add_genome_to_config,
                "download_genome_files": tool_module.download_genome_files,
            }
            
            # --- 修正后的逻辑：优先处理特殊工具 ---
            if function_name == "get_task_status":
                task_id = function_args.get("task_id")
                with db_lock:
                    task_info = TASK_DATABASE.get(task_id)
                
                if not task_info:
                    function_response = f"错误：找不到任务 '{task_id}'。"
                else:
                    pid = task_info.get("pid")
                    current_status = task_info.get("status", "unknown")
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
                                    TASK_DATABASE[task_id]["status"] = "finished"
                                    TASK_DATABASE[task_id]["exit_code"] = exit_code
                        except OSError:
                            process_status_text = "已结束 (进程不存在)"
                            with db_lock:
                                if TASK_DATABASE.get(task_id):
                                    TASK_DATABASE[task_id]["status"] = "finished"
                    elif current_status == "finished":
                        exit_code = task_info.get('exit_code', 'N/A')
                        process_status_text = f"已结束 (退出码: {exit_code})"
                    function_response = f"任务 '{task_id}' 的状态为: {process_status_text} (PID: {pid})。"
            
            # --- 然后处理通用工具 ---
            else:
                function_to_call = available_tools.get(function_name)
                if not function_to_call:
                    function_response = f"错误：未知的工具名称 '{function_name}'"
                else:
                    try:
                        # 统一调用工具
                        tool_result = function_to_call(**function_args)
                        
                        # 如果是运行流程的工具，需要特殊处理结果以创建任务
                        if function_name == "run_rna_seq_pipeline":
                            if tool_result.get("status") == "success":
                                global task_id_counter
                                with db_lock:
                                    task_id = f"task_{task_id_counter}"
                                    task_id_counter += 1
                                    TASK_DATABASE[task_id] = {
                                        "pid": tool_result.get("pid"),
                                        "srr_list_file": tool_result.get("srr_list_file"),
                                        "status": "running",
                                        "start_time": time.time(),
                                        "srr_list": function_args.get("srr_list")
                                    }
                                function_response = f"任务 '{task_id}' 已成功启动。进程 PID: {tool_result.get('pid')}。"
                            else:
                                function_response = f"启动流程失败: {tool_result.get('message', '未知工具错误')}"
                        else:
                            # 对于其他所有工具，直接将字典结果转为字符串
                            function_response = json.dumps(tool_result, ensure_ascii=False, indent=4)

                    except Exception as e:
                        function_response = f"执行工具 '{function_name}' 时发生错误: {e}"

            # 5. 将工具执行结果添加到消息历史中
            messages.append(
                {
                    "tool_call_id": tool_call.id,
                    "role": "tool",
                    "name": function_name,
                    "content": function_response,
                }
            )
        
        # 6. 第二次调用 LLM，让它根据工具结果生成最终回复
        print(f"--- 第二次调用 LLM (携带工具结果) ---")
        print(f"发送的消息: {messages}")
        
        try:
            second_response = client.chat.completions.create(
                model=model_name,
                messages=messages,
                stream=True, # 以流式模式获取最终回复
            )
            # 7. 流式传输最终回复
            for chunk in second_response:
                content = chunk.choices[0].delta.content
                if content:
                    response_chunk = {
                        "id": chunk.id, "object": chunk.object, "created": chunk.created, "model": chunk.model,
                        "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": None}]
                    }
                    yield f"data: {json.dumps(response_chunk)}\n\n"
        except Exception as e:
            print(f"第二次调用 LLM API 时发生错误: {e}")
            yield f"data: {json.dumps({'error': str(e)})}\n\n"

    else:
        # 如果 LLM 不需要调用工具，直接流式返回它的回复
        print("--- LLM 无需调用工具，直接回复 ---")
        content = response_message.content or ""
        for char in content:
            response_chunk = {
                "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                "choices": [{"index": 0, "delta": {"content": char}, "finish_reason": None}]
            }
            yield f"data: {json.dumps(response_chunk)}\n\n"
            await asyncio.sleep(0.02)

    # 发送结束标志
    final_chunk = {
        "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
        "choices": [{"index": 0, "delta": {}, "finish_reason": "stop"}]
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