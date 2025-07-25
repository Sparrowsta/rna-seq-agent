# agent/server.py

import json
import time
import asyncio
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import List, AsyncGenerator
from fastapi.middleware.cors import CORSMiddleware

# --- 1. 定义数据模型 ---
class StandardMessage(BaseModel):
    role: str
    content: str

class ChatInput(BaseModel):
    messages: List[StandardMessage]

# --- 2. 创建 FastAPI 应用实例 ---
app = FastAPI(
    title="从零构建的流式 Agent 服务器",
    description="这是一个支持流式响应、符合OpenAI接口标准的Agent服务器。",
    version="0.2.4", # 版本更新
)

# --- 3. 添加 CORS 中间件 ---
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

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

async def stream_generator() -> AsyncGenerator[str, None]:
    """
    异步生成器，用于模拟流式响应。
    """
    model_id = "rna-seq-agent-v1"
    response_chunks = [
        "Agent ", "已收到", "您的", "消息。", " ",
        "我", "现在", "是一个", "流式", "Agent。",
    ]
    
    for chunk in response_chunks:
        response_chunk = {
            "id": f"chatcmpl-{int(time.time())}",
            "object": "chat.completion.chunk",
            "created": int(time.time()),
            "model": model_id,
            "choices": [{
                "index": 0,
                "delta": {"content": chunk},
                "finish_reason": None
            }]
        }
        sse_formatted_chunk = f"data: {json.dumps(response_chunk)}\n\n"
        yield sse_formatted_chunk
        await asyncio.sleep(0.1)

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
    核心聊天接口，返回一个流式响应。
    """
    return StreamingResponse(stream_generator(), media_type="text/event-stream")

# --- 5. 启动服务器 ---
if __name__ == "__main__":
    import uvicorn
    print("正在启动流式 Agent 服务器 v0.2.4 (移除/v1前缀)，访问 http://localhost:8001")
    uvicorn.run(app, host="0.0.0.0", port=8001)