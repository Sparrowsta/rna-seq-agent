# agent/server.py

import json
import time
import asyncio
import os
import tempfile
import re
import json
from enum import Enum
import hashlib
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
# --- 会话状态管理 ---
class SessionState(str, Enum):
    CONVERSING = "CONVERSING"
    ANALYZING = "ANALYZING"

# 使用一个全局字典来模拟会话存储
# key 是会话ID, value 是状态
session_states: Dict[str, SessionState] = {}

import agent.tools as tool_module # 导入整个模块以便于函数查找

# --- 1. 定义数据模型 ---
class StandardMessage(BaseModel):
    role: str
    content: str

class ChatInput(BaseModel):
    messages: List[StandardMessage]


def get_session_id(chat_input: ChatInput) -> str:
    """
    根据聊天记录生成一个稳定的会话ID。
    """
    # 使用所有消息内容的哈希值作为ID
    # 注意：这只是一个简单的实现，在生产环境中需要更健壮的会话管理
    message_content = "".join([msg.content for msg in chat_input.messages])
    return f"session_{hashlib.md5(message_content.encode()).hexdigest()}"


# --- 2. 创建 FastAPI 应用实例 ---
app = FastAPI(
    title="LLM 驱动的流式 Agent 服务器",
    description="这是一个由大型语言模型驱动、支持工具调用、符合OpenAI接口标准的Agent服务器。",
    version="1.0.0", # 版本跃迁
)

# --- 3. 配置和初始化 ---
# 加载 .env 文件中的环境变量
# 这应该在访问任何环境变量之前完成
# 加载项目根目录的 .env 文件 (如果存在)，
# 这使得在容器外本地运行也成为可能。
# 在容器内，这些变量主要由 docker-compose 的 env_file 指令注入。
load_dotenv()

# 初始化 OpenAI 客户端
# 从环境变量中获取配置
api_key = os.getenv("OPENAI_API_KEY")
base_url = os.getenv("OPENAI_API_BASE")
model_name = os.getenv("OPENAI_MODEL_NAME", "default-model") 

if not api_key or not base_url:
    raise ValueError("请在项目根目录的 .env 文件中设置 OPENAI_API_KEY 和 OPENAI_API_BASE")

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
    一个由 LLM 驱动的、支持状态机和工具调用的 Agent 响应生成器。
    """
    session_id = get_session_id(chat_input)
    
    # 1. 获取当前会话状态，默认为 CONVERSING
    current_state = session_states.get(session_id, SessionState.CONVERSING)
    print(f"--- 会话 {session_id}: 当前状态 {current_state.value} ---")

    # 2. 准备发送给 LLM 的消息
    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    for msg in chat_input.messages:
        if msg.role != "system":
            messages.append({"role": msg.role, "content": msg.content})

    # 3. 动态注入当前状态到最新的用户消息中
    state_message = f"[session_state: {current_state.value}]"
    if messages[-1]["role"] == "user":
        messages[-1]["content"] = f"{messages[-1]['content']}\n\n{state_message}"
    else:
        # 如果最后一条消息不是用户消息，则添加一条系统消息来传递状态
        messages.append({"role": "system", "content": state_message})

    # 4. 根据状态决定执行路径
    while True: # 使用循环来处理状态转换
        if current_state == SessionState.CONVERSING:
            # --- 对话模式 ---
            print(f"--- 会话 {session_id}: 进入 CONVERSING 模式 ---")
            
            # 直接调用LLM进行一次对话，检查是否要切换状态
            response = client.chat.completions.create(
                model=model_name,
                messages=messages,
                tools=TOOLS,
                tool_choice="auto", # 允许模型自主决定是否调用工具
                temperature=0.0,
                stream=False  # 在对话模式下，我们需要先获得完整响应来判断是否切换状态
            )
            
            response_message = response.choices[0].message
            tool_calls = response_message.tool_calls

            # 检查LLM是否调用了 'start_analysis_tool'
            if tool_calls and any(tc.function.name == 'start_analysis_tool' for tc in tool_calls):
                session_states[session_id] = SessionState.ANALYZING
                current_state = SessionState.ANALYZING # 更新本地状态
                print(f"--- 会话 {session_id}: 状态切换 -> ANALYZING ---")
                
                # 向前端发送状态切换的确认信息
                switch_message = "✅ 好的，已确认分析计划。状态已切换到分析模式，我将开始执行任务。请注意，在分析完成前我将专注于执行，无法进行新的对话。\n\n"
                switch_chunk = {
                    "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                    "choices": [{"index": 0, "delta": {"content": switch_message}, "finish_reason": None}] # finish_reason is None because we continue
                }
                yield f"data: {json.dumps(switch_chunk)}\n\n"
                
                # 不要返回，继续循环以进入分析模式
                continue

            else:
                # 正常流式返回对话内容
                # 由于上面已经进行了一次非流式调用，我们模拟流式返回
                content = response_message.content or ""
                if content:
                    response_chunk = {
                        "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                        "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": None}]
                    }
                    yield f"data: {json.dumps(response_chunk)}\n\n"
                
                # 如果有其他工具调用（非start_analysis_tool），也一并处理
                if tool_calls:
                    messages.append(response_message)
                    
                    # --- 执行工具调用 ---
                    available_tools = {
                        "search_genome_tool": tool_module.search_genome_tool,
                        "search_fastq_tool": tool_module.search_fastq_tool,
                        "check_files_exist_tool": tool_module.check_files_exist_tool,
                        "list_files": tool_module.list_files,
                        "get_task_status": tool_module.get_task_status,
                        "check_environment_tool": tool_module.check_environment_tool,
                    }

                    for tool_call in tool_calls:
                        function_name = tool_call.function.name
                        function_args = json.loads(tool_call.function.arguments)
                        function_response = ""
                        
                        try:
                            function_to_call = available_tools.get(function_name)
                            if not function_to_call:
                                function_response = f"错误: 在对话模式下找不到或不允许使用工具 '{function_name}'。"
                            else:
                                # 准备工具函数参数，处理特殊依赖
                                tool_kwargs = function_args.copy()
                                if function_name == "get_task_status":
                                    tool_kwargs["task_database"] = TASK_DATABASE
                                    tool_kwargs["db_lock"] = db_lock
                                
                                # 调用工具函数
                                print(f"Calling tool: {function_name} with args: {tool_kwargs}")
                                tool_result = function_to_call(**tool_kwargs)
                                print(f"Tool {function_name} returned: {tool_result}")
                                
                                # 将结果格式化为字符串
                                if isinstance(tool_result, dict):
                                    function_response = json.dumps(tool_result, ensure_ascii=False)
                                else:
                                    function_response = str(tool_result)

                        except Exception as e:
                            function_response = f"执行工具 '{function_name}' 时出错: {e}"

                        messages.append({
                            "tool_call_id": tool_call.id,
                            "role": "tool",
                            "name": function_name,
                            "content": function_response,
                        })

                    # --- 工具调用后，再次调用LLM以获得最终的自然语言响应 ---
                    final_response_stream = client.chat.completions.create(
                        model=model_name,
                        messages=messages,
                        temperature=0.0,
                        stream=True
                    )

                    # 将最终响应流式传输回客户端
                    for chunk in final_response_stream:
                        chunk_data = chunk.model_dump_json()
                        yield f"data: {chunk_data}\n\n"
                    
                # 对话或普通工具调用完成后，结束流程
                break

        elif current_state == SessionState.ANALYZING:
            # --- 分析模式 (React循环) ---
            print(f"--- 会话 {session_id}: 进入 ANALYZING 模式 ---")
            react_cycle_count = 0
            max_react_cycles = 100

            while react_cycle_count < max_react_cycles:
                react_cycle_count += 1
                
                # 调用 LLM
                try:
                    response = client.chat.completions.create(
                        model=model_name,
                        messages=messages,
                        tools=TOOLS,
                        temperature=0.0
                    )
                except Exception as e:
                    yield f"data: {json.dumps({'error': str(e)})}\n\n"
                    break

                response_message = response.choices[0].message
                tool_calls = response_message.tool_calls

                # 将LLM的思考或行动请求加入历史
                messages.append(response_message)

                if tool_calls:
                    # 执行工具调用
                    for tool_call in tool_calls:
                        function_name = tool_call.function.name
                        function_args = json.loads(tool_call.function.arguments)
                        
                        # --- 详细的工具调用逻辑 ---
                        available_tools = {
                            "get_task_status": tool_module.get_task_status,
                            "list_files": tool_module.list_files,
                            "add_genome_to_config": tool_module.add_genome_to_config,
                            "unsupported_request": tool_module.unsupported_request,
                            "check_environment_tool": tool_module.check_environment_tool,
                            "setup_environment_tool": tool_module.setup_environment_tool,
                            "search_genome_tool": tool_module.search_genome_tool,
                            "search_fastq_tool": tool_module.search_fastq_tool,
                            "download_genome_tool": tool_module.download_genome_tool,
                            "download_fastq_tool": tool_module.download_fastq_tool,
                            "validate_fastq_tool": tool_module.validate_fastq_tool,
                            "check_files_exist_tool": tool_module.check_files_exist_tool,
                            "build_star_index_tool": tool_module.build_star_index_tool,
                            "run_fastp_tool": tool_module.run_fastp_tool,
                            "run_star_align_tool": tool_module.run_star_align_tool,
                            "run_featurecounts_tool": tool_module.run_featurecounts_tool,
                            "collect_results_tool": tool_module.collect_results_tool,
                            "generate_report_tool": tool_module.generate_report_tool,
                            "start_analysis_tool": tool_module.start_analysis_tool,
                            "react_status_tool": tool_module.react_status_tool,
                            "react_plan_tool": tool_module.react_plan_tool,
                            "react_evaluate_tool": tool_module.react_evaluate_tool,
                            "react_summary_tool": tool_module.react_summary_tool,
                            "validate_tool_call_format": tool_module.validate_tool_call_format,
                        }
                        
                        global task_id_counter
                        function_response = ""

                        try:
                            function_to_call = available_tools.get(function_name)
                            if not function_to_call:
                                function_response = f"错误：未知的工具名称 '{function_name}'"
                            else:
                                tool_kwargs = function_args.copy()
                                TOOLS_NEEDING_DB = {"get_task_status"}
                                TOOLS_NEEDING_COUNTER = {}

                                if function_name in TOOLS_NEEDING_DB:
                                    tool_kwargs["task_database"] = TASK_DATABASE
                                    tool_kwargs["db_lock"] = db_lock
                                
                                if function_name in TOOLS_NEEDING_COUNTER:
                                    tool_kwargs["task_id_counter"] = task_id_counter

                                print(f"Calling tool: {function_name} with args: {tool_kwargs}")
                                tool_result = function_to_call(**tool_kwargs)
                                print(f"Tool {function_name} returned: {tool_result}")

                                if function_name in TOOLS_NEEDING_COUNTER and isinstance(tool_result, dict):
                                    with db_lock:
                                        task_id_counter = tool_result.get("updated_task_id_counter", task_id_counter)

                                if isinstance(tool_result, dict):
                                    function_response = json.dumps(tool_result, ensure_ascii=False, indent=2)
                                    display_message = tool_result.get("message", f"状态: {tool_result.get('status', '执行完成')}")
                                else:
                                    function_response = str(tool_result)
                                    display_message = function_response

                                tool_result_chunk = {
                                    "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                                    "choices": [{"index": 0, "delta": {"content": f"✅ 工具 {function_name} 执行完成: {display_message}\n\n"}, "finish_reason": None}]
                                }
                                yield f"data: {json.dumps(tool_result_chunk)}\n\n"

                        except Exception as e:
                            function_response = f"执行工具 '{function_name}' 时发生严重错误: {e}"
                            tool_error_chunk = {
                                "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                                "choices": [{"index": 0, "delta": {"content": f"❌ 工具 {function_name} 执行失败: {e}"}, "finish_reason": None}]
                            }
                            yield f"data: {json.dumps(tool_error_chunk)}\n\n"

                        messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})

                        # 检查任务是否完成或失败
                        task_finished = (function_name == 'generate_report_tool' and 'error' not in function_response.lower())
                        task_failed = 'error' in function_response.lower()

                        if task_finished or task_failed:
                            session_states[session_id] = SessionState.CONVERSING
                            print(f"--- 会话 {session_id}: 分析完成/失败，状态切换 -> CONVERSING ---")
                            
                            final_message = "分析流程已结束。您可以提出新问题或开始新的分析任务。"
                            if task_failed:
                                final_message = f"分析流程因错误而终止。错误信息: {function_response}。请检查问题后重试。"

                            final_chunk = {
                                "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                                "choices": [{"index": 0, "delta": {"content": f"\n{final_message}"}, "finish_reason": "stop"}]
                            }
                            yield f"data: {json.dumps(final_chunk)}\n\n"
                            yield "data: [DONE]\n\n"
                            return # 明确返回以终止生成器

                else: # LLM没有调用工具，只是返回文本
                    content = response_message.content or ""
                    if content:
                        response_chunk = {
                            "id": response.id, "object": "chat.completion.chunk", "created": response.created, "model": response.model,
                            "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": None}]
                        }
                        yield f"data: {json.dumps(response_chunk)}\n\n"
                    
                    # 如果在分析模式下LLM只是说话，通常意味着分析结束了
                    session_states[session_id] = SessionState.CONVERSING
                    print(f"--- 会话 {session_id}: 分析完成，状态切换 -> CONVERSING ---")
                    break # 退出内部 react 循环
            
            # 分析流程已通过自然语言响应结束，现在退出主循环
            break

    # 发送结束标志
    final_chunk = {
        "id": "final_chunk", "object": "chat.completion.chunk", "created": int(time.time()), "model": model_name,
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




@app.get("/task-status/{task_id}")
def get_task_status_endpoint(task_id: str):
    """用于轮询长时任务状态的端点。"""
    print(f"--- [状态端点] 正在查询任务: {task_id} ---")
    try:
        # 直接调用我们的工具函数，并传入所需的依赖
        status_result = tool_module.get_task_status(
            task_id=task_id,
            task_database=TASK_DATABASE,
            db_lock=db_lock
        )
        # 如果工具函数返回了它自己的错误（例如 "找不到任务"），则将其转换为 404
        if not status_result or status_result.get("status") == "error":
            raise HTTPException(status_code=404, detail=status_result.get("message", f"找不到任务 '{task_id}'。"))
        
        return status_result
    except HTTPException as he:
        raise he # 重新抛出已有的 HTTP 异常
    except Exception as e:
        print(f"!!! [状态端点] 查询状态时出错: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# --- 5. 启动服务器 ---
if __name__ == "__main__":
    import uvicorn
    print("正在启动流式 Agent 服务器 v0.3.0 ，访问 http://localhost:48001")
    uvicorn.run(app, host="0.0.0.0", port=48001)