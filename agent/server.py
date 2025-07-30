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
    一个由 LLM 驱动的、支持工具调用的 Agent 响应生成器。
    现在支持React模式：思考-行动-观察循环。
    """
    # 检查是否为首次交互（只有一条用户消息且没有历史记录）
    is_first_interaction = len(chat_input.messages) == 1 and chat_input.messages[0].role == "user"
    
    # 1. 准备发送给 LLM 的消息，确保我们的系统提示是唯一的
    messages = [{"role": "system", "content": SYSTEM_PROMPT}]
    # 过滤掉任何可能从上游传入的 system 消息，只保留 user 和 assistant 的消息
    for msg in chat_input.messages:
        if msg.role != "system":
            messages.append({"role": msg.role, "content": msg.content})
    
    # 如果是首次交互，发送欢迎信息
    if is_first_interaction:
        welcome_message = """👋 欢迎使用RNA-seq 分析平台

您可以这样开始：
1. 输入SRR号和参考基因组名，例如：`帮我分析SRR17469059 基因组用mm10`
2. 系统会自动完成数据下载、质量控制、比对、定量和报告生成。
3. 支持智能跳过已完成步骤，节省计算资源。

常用命令：
- 查看可用基因组：`列出可用基因组`
- 添加新基因组：`添加基因组 mm10 物种 mouse ...`
- 查询分析进度：`查询任务状态`

祝您分析顺利！"""
        
        # 发送欢迎信息
        welcome_chunk = {
            "id": "welcome",
            "object": "chat.completion.chunk",
            "created": int(time.time()),
            "model": "welcome",
            "choices": [{"index": 0, "delta": {"content": welcome_message}, "finish_reason": None}]
        }
        yield f"data: {json.dumps(welcome_chunk)}\n\n"
        
        # 添加欢迎信息到消息历史
        messages.append({"role": "assistant", "content": welcome_message})

    # React模式状态跟踪
    react_cycle_count = 0
    max_react_cycles = 100  # 防止无限循环
    
    while react_cycle_count < max_react_cycles:
        react_cycle_count += 1
        
        # 2. 调用 LLM，让它进行思考并决定行动
        print(f"--- React循环 {react_cycle_count}: 调用 LLM (模型: {model_name}) ---")
        print(f"发送的消息: {messages}")
        print(f"可用的工具: {TOOLS}")
        
        try:
            response = client.chat.completions.create(
                model=model_name,
                messages=messages,
                tools=TOOLS,
                temperature=0.0
            )
        except Exception as e:
            print(f"调用 LLM API 时发生错误: {e}")
            yield f"data: {json.dumps({'error': str(e)})}\n\n"
            yield "data: [DONE]\n\n"
            return

        response_message = response.choices[0].message
        tool_calls = response_message.tool_calls

            # 3. 发送思考阶段的信息
        thought_chunk = {
            "id": response.id, 
            "object": "chat.completion.chunk", 
            "created": response.created, 
            "model": response.model,
            "choices": [{"index": 0, "delta": {"content": f"🤔 思考阶段 (循环 {react_cycle_count}): LLM正在分析当前情况并制定行动计划...\n"}, "finish_reason": None}]
            }
        yield f"data: {json.dumps(thought_chunk)}\n\n"
            
            # 强制立即显示
        yield f"data: {json.dumps({'id': response.id, 'object': 'chat.completion.chunk', 'created': response.created, 'model': response.model, 'choices': [{'index': 0, 'delta': {'content': ''}, 'finish_reason': None}]})}\n\n"

        # 4. 检查 LLM 是否决定调用工具
        if tool_calls:
            print(f"--- React循环 {react_cycle_count}: LLM 决定调用工具: {tool_calls} ---")
            messages.append(response_message)  # 将 assistant 的回复（包括工具调用请求）添加到历史记录中

            # 5. 执行所有工具调用（行动阶段）
            action_chunk = {
                "id": response.id, 
                "object": "chat.completion.chunk", 
                "created": response.created, 
                "model": response.model,
                "choices": [{"index": 0, "delta": {"content": f"🔧 行动阶段 (循环 {react_cycle_count}): 执行 {len(tool_calls)} 个工具调用..."}, "finish_reason": None}]
            }
            yield f"data: {json.dumps(action_chunk)}\n\n"

            for tool_call in tool_calls:
                function_name = tool_call.function.name
                
                try:
                    function_args = json.loads(tool_call.function.arguments)
                except json.JSONDecodeError:
                    print(f"❌ JSON解析错误: {tool_call.function.arguments}")
                    function_response = f"错误: LLM 返回了无效的 JSON 参数: {tool_call.function.arguments}"
                    messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                    continue
                
                # 检查是否有特殊token或错误格式
                if "REDACTED_SPECIAL_TOKEN" in tool_call.function.arguments:
                    print(f"❌ 检测到特殊token: {tool_call.function.arguments}")
                    function_response = f"错误: 检测到特殊token，请使用标准的工具调用格式"
                    messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                    continue
                
                # 检查是否有function标签
                if "function" in tool_call.function.arguments:
                    print(f"❌ 检测到function标签: {tool_call.function.arguments}")
                    function_response = f"错误: 检测到function标签，请使用标准的工具调用格式"
                    messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                    continue
                
                # 检查是否有JSON标签
                if "<JSON>" in tool_call.function.arguments:
                    print(f"❌ 检测到JSON标签: {tool_call.function.arguments}")
                    function_response = f"错误: 检测到JSON标签，请使用标准的工具调用格式"
                    messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                    continue
                
                # 尝试清理和修复参数格式
                cleaned_args = tool_call.function.arguments
                if "REDACTED_SPECIAL_TOKEN" in cleaned_args:
                    # 尝试提取JSON部分
                    import re
                    json_match = re.search(r'\{[^}]*\}', cleaned_args)
                    if json_match:
                        cleaned_args = json_match.group(0)
                        print(f"🔧 尝试清理参数: {cleaned_args}")
                        try:
                            function_args = json.loads(cleaned_args)
                        except json.JSONDecodeError:
                            print(f"❌ 清理后的参数仍然无效: {cleaned_args}")
                            function_response = f"错误: 无法解析工具参数，请重新调用工具"
                            messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                            continue
                    else:
                        print(f"❌ 无法从参数中提取有效JSON: {cleaned_args}")
                        function_response = f"错误: 无法解析工具参数，请重新调用工具"
                        messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                        continue
                
                # 检查是否包含多个REDACTED_SPECIAL_TOKEN（更严格的检查）
                if cleaned_args.count("REDACTED_SPECIAL_TOKEN") > 1:
                    print(f"❌ 检测到多个REDACTED_SPECIAL_TOKEN: {cleaned_args}")
                    function_response = f"错误: 检测到多个REDACTED_SPECIAL_TOKEN，请使用标准的工具调用格式"
                    messages.append({"tool_call_id": tool_call.id, "role": "tool", "name": function_name, "content": function_response})
                    continue

                print(f"🔧 执行工具: {function_name}")
                print(f"📝 参数: {function_args}")
                print(f"🔍 参数类型: {type(function_args)}")

                # 实时输出工具调用信息
                tool_call_chunk = {
                    "id": response.id, 
                    "object": "chat.completion.chunk", 
                    "created": response.created, 
                    "model": response.model,
                    "choices": [{"index": 0, "delta": {"content": f"🔧 调用工具: {function_name} (参数: {function_args})\n"}, "finish_reason": None}]
                }
                yield f"data: {json.dumps(tool_call_chunk)}\n\n"
                
                # 强制立即显示
                yield f"data: {json.dumps({'id': response.id, 'object': 'chat.completion.chunk', 'created': response.created, 'model': response.model, 'choices': [{'index': 0, 'delta': {'content': ''}, 'finish_reason': None}]})}\n\n"

                # --- 最终重构：完全统一的工具调用逻辑 ---
                # --- 新的、模块化的工具集 ---
                available_tools = {
                    # v5.2 Tools
                    "get_task_status": tool_module.get_task_status,
                    "list_files": tool_module.list_files,
                    "add_genome_to_config": tool_module.add_genome_to_config,
                    "unsupported_request": tool_module.unsupported_request,
                    # React模式工具
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
                        # --- 新的、更清晰的依赖注入逻辑 ---
                        tool_kwargs = function_args.copy()

                        # 定义哪些工具需要哪些共享资源
                        TOOLS_NEEDING_DB = {
                            "get_task_status",
                        }
                        TOOLS_NEEDING_COUNTER = {}

                        # 注入数据库和锁
                        if function_name in TOOLS_NEEDING_DB:
                            tool_kwargs["task_database"] = TASK_DATABASE
                            tool_kwargs["db_lock"] = db_lock
                        
                        # 注入任务ID计数器
                        if function_name in TOOLS_NEEDING_COUNTER:
                            # 传递当前的计数器值
                            tool_kwargs["task_id_counter"] = task_id_counter

                        # 调用工具
                        tool_result = function_to_call(**tool_kwargs)

                        # 如果工具更新了计数器，则同步回全局计数器
                        # 这是为了让 create_task 和 download_genome 都能安全地增加ID
                        if function_name in TOOLS_NEEDING_COUNTER and isinstance(tool_result, dict):
                            with db_lock:
                                task_id_counter = tool_result.get("updated_task_id_counter", task_id_counter)

                        # 统一处理返回结果
                        # 如果结果是字典，序列化为 JSON 字符串
                        if isinstance(tool_result, dict):
                            function_response = json.dumps(tool_result, ensure_ascii=False, indent=2)
                            # 提取关键信息用于显示
                            if "message" in tool_result:
                                display_message = tool_result["message"]
                            elif "status" in tool_result:
                                display_message = f"状态: {tool_result['status']}"
                            else:
                                display_message = "执行完成"
                        else:
                            # 否则，直接使用返回的字符串
                            function_response = str(tool_result)
                            display_message = function_response

                        # 实时输出工具执行结果
                        tool_result_chunk = {
                            "id": response.id, 
                            "object": "chat.completion.chunk", 
                            "created": response.created, 
                            "model": response.model,
                            "choices": [{"index": 0, "delta": {"content": f"✅ 工具 {function_name} 执行完成: {display_message}\n\n"}, "finish_reason": None}]
                        }
                        yield f"data: {json.dumps(tool_result_chunk)}\n\n"
                        
                        # 强制刷新输出
                        yield f"data: {json.dumps({'id': response.id, 'object': 'chat.completion.chunk', 'created': response.created, 'model': response.model, 'choices': [{'index': 0, 'delta': {'content': ''}, 'finish_reason': None}]})}\n\n"

                except Exception as e:
                    # 统一的异常处理
                    function_response = f"执行工具 '{function_name}' 时发生严重错误: {e}"
                    
                    # 实时输出工具执行错误
                    tool_error_chunk = {
                        "id": response.id, 
                        "object": "chat.completion.chunk", 
                        "created": response.created, 
                        "model": response.model,
                        "choices": [{"index": 0, "delta": {"content": f"❌ 工具 {function_name} 执行失败: {e}"}, "finish_reason": None}]
                    }
                    yield f"data: {json.dumps(tool_error_chunk)}\n\n"

                # 6. 将工具执行结果添加到消息历史中
                messages.append(
                    {
                        "tool_call_id": tool_call.id,
                        "role": "tool",
                        "name": function_name,
                        "content": function_response,
                    }
                )
            
            # 7. 观察阶段 - 分析工具执行结果
            observation_chunk = {
                "id": response.id, 
                "object": "chat.completion.chunk", 
                "created": response.created, 
                "model": response.model,
                "choices": [{"index": 0, "delta": {"content": f"👀 观察阶段 (循环 {react_cycle_count}): 正在分析工具执行结果并评估下一步行动...\n"}, "finish_reason": None}]
            }
            yield f"data: {json.dumps(observation_chunk)}\n\n"
            
            # 强制立即显示
            yield f"data: {json.dumps({'id': response.id, 'object': 'chat.completion.chunk', 'created': response.created, 'model': response.model, 'choices': [{'index': 0, 'delta': {'content': ''}, 'finish_reason': None}]})}\n\n"
            
            # 8. 检查是否需要继续React循环
            # 分析工具执行结果，决定是否继续循环
            should_continue = False
            
            for tool_call in tool_calls:
                function_name = tool_call.function.name
                
                # 检查是否是计划工具，如果是，需要继续执行
                if function_name in ["plan_analysis_task", "react_plan_tool"]:
                    should_continue = True
                    break
                
                # 检查是否有需要进一步处理的工具（如长时间运行的任务）
                if function_name in ["download_genome_files"]:
                    should_continue = True
                    break
                
                # 检查是否是查询工具，如果是，需要继续执行（让LLM决定下一步）
                if function_name in ["get_task_status", "list_available_genomes", "list_files", "search_genome_tool", "search_fastq_tool", "validate_fastq_tool", "check_environment_tool", "check_files_exist_tool", "build_star_index_tool"]:
                    should_continue = True
                    break
                
                # 检查是否是分析工具，如果是，需要继续执行（让LLM决定下一步）
                if function_name in ["run_fastp_tool", "run_star_align_tool", "run_featurecounts_tool", "collect_results_tool"]:
                    should_continue = True
                    break
                
                # 检查是否是最终工具，如果是，可以结束循环
                if function_name in ["generate_report_tool", "react_summary_tool"]:
                    should_continue = False
                    break
            
            if should_continue:
                # 继续React循环，让LLM决定下一步行动
                print(f"--- React循环 {react_cycle_count}: 继续循环，等待LLM决定下一步 ---")
                continue
            else:
                # 所有工具都执行完成，准备生成最终回复
                print(f"--- React循环 {react_cycle_count}: 所有工具执行完成，准备生成最终回复 ---")
                break

        else:
            # LLM没有调用工具，可能只是输出思考
            print(f"--- React循环 {react_cycle_count}: LLM 输出内容，无工具调用 ---")
            content = response_message.content or ""
            
            # 将LLM的思考内容以OpenAI格式发送到前端，但不终止React循环
            if content.strip():
                # 以OpenAI格式发送思考内容，让CherryStudio能够正常显示
                response_chunk = {
                    "id": response.id, 
                    "object": "chat.completion.chunk", 
                    "created": response.created, 
                    "model": response.model,
                    "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": None}]
                }
                yield f"data: {json.dumps(response_chunk)}\n\n"
            
            # 将LLM的回复添加到消息历史中
            messages.append(response_message)
            
            # 检查内容是否包含明确的结束信号
            if any(keyword in content.lower() for keyword in ["分析完成", "流程结束", "任务完成", "报告生成完成", "最终总结", "分析报告"]):
                print(f"--- React循环 {react_cycle_count}: 检测到完成信号，结束循环 ---")
                break
            else:
                # 继续React循环，让LLM决定下一步行动
                print(f"--- React循环 {react_cycle_count}: 继续循环，等待LLM决定下一步 ---")
                continue

    # 9. 最终回复阶段 - 生成总结性回复
    if react_cycle_count > 0:
        print(f"--- 生成最终回复 (React循环完成: {react_cycle_count}) ---")
        print(f"发送的消息: {messages}")
        
        try:
            final_response = client.chat.completions.create(
                model=model_name,
                messages=messages,
                stream=True, # 以流式模式获取最终回复
                temperature=0.0
            )
            # 流式传输最终回复
            for chunk in final_response:
                content = chunk.choices[0].delta.content
                if content:
                    response_chunk = {
                        "id": chunk.id, "object": chunk.object, "created": chunk.created, "model": chunk.model,
                        "choices": [{"index": 0, "delta": {"content": content}, "finish_reason": None}]
                    }
                    yield f"data: {json.dumps(response_chunk)}\n\n"
        except Exception as e:
            print(f"生成最终回复时发生错误: {e}")
            yield f"data: {json.dumps({'error': str(e)})}\n\n"

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