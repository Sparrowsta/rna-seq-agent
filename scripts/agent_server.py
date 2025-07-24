import os
import json
from fastapi import FastAPI
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import List, Optional

from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langchain.agents import AgentExecutor, create_react_agent
from langchain_core.prompts import ChatPromptTemplate
from scripts.agent import tools, prompt_template
from mcpserver.server import app as mcpserver_app

# The .env file is mounted at /app/config/.env
dotenv_path = '/app/config/.env'
if os.path.exists(dotenv_path):
    load_dotenv(dotenv_path=dotenv_path)
    print(f"Loaded .env file from {dotenv_path}")
else:
    print(f"Warning: .env file not found at {dotenv_path}")

app = FastAPI(
    title="Agent Server",
    description="An OpenAI-compatible server for the bioinformatics agent, which also hosts the mcpserver.",
    version="1.0.0"
)

# Mount the mcpserver app under the /mcpserver path
app.mount("/mcpserver", mcpserver_app)


# --- Agent Initialization ---
llm = ChatOpenAI(
    temperature=0,
    model_name=os.environ.get("OPENAI_MODEL_NAME"),
    openai_api_key=os.environ.get("OPENAI_API_KEY"),
    openai_api_base=os.environ.get("OPENAI_API_BASE"),
)
prompt = ChatPromptTemplate.from_template(prompt_template)
agent = create_react_agent(llm, tools, prompt)
agent_executor = AgentExecutor(agent=agent, tools=tools, verbose=True, handle_parsing_errors=True)

# --- Pydantic Models for OpenAI-like requests ---
class StandardMessage(BaseModel):
    role: str
    content: str

class ChatInput(BaseModel):
    messages: List[StandardMessage]
    # We can ignore other fields like model, temperature, stream for now

@app.post("/v1/chat/completions")
async def chat_completions(chat_input: ChatInput):
    """
    Receives a message from the user and returns the agent's response as a stream.
    The path is now /v1/chat/completions to be more compliant.
    """
    last_user_message_text = ""
    # Iterate backwards to find the last user message
    for message in reversed(chat_input.messages):
        if message.role == 'user':
            last_user_message_text = message.content
            break
    
    if not last_user_message_text:
        last_user_message_text = ""

    async def event_generator():
        full_response = ""
        # Use astream for async iteration
        async for chunk in agent_executor.astream({"input": last_user_message_text}):
            # We are interested in the 'output' chunks from the final answer
            if "output" in chunk:
                content = chunk["output"]
                # The output may come in chunks or as a whole.
                # We compare with the previously sent content to find the new part.
                new_content = content.replace(full_response, "", 1)
                if new_content:
                    full_response = content
                    response_chunk = {
                        "choices": [{
                            "delta": {"content": new_content},
                            "index": 0,
                        }]
                    }
                    yield f"data: {json.dumps(response_chunk)}\n\n"
        
        # Send the final empty chunk with finish_reason
        final_chunk = {
            "choices": [{
                "delta": {},
                "index": 0,
                "finish_reason": "stop"
            }]
        }
        yield f"data: {json.dumps(final_chunk)}\n\n"
        yield "data: [DONE]\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")

@app.get("/v1/models")
async def list_models():
    """
    Returns a list of available models, mimicking the OpenAI models endpoint.
    The path is now /v1/models to be more compliant.
    """
    model_name = os.environ.get("OPENAI_MODEL_NAME", "default-model")
    return {
        "object": "list",
        "data": [
            {
                "id": model_name,
                "object": "model",
                "created": 1677610602,
                "owned_by": "system",
            }
        ],
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)