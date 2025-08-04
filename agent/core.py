import os
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from .tools import list_directory

tools = [list_directory]
llm = ChatGoogleGenerativeAI(
    model=os.environ.get("GEMINI_MODEL_NAME"),
    google_api_key=os.environ.get("GOOGLE_API_KEY")
)
llm_with_tools = llm.bind_tools(tools)

prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            "You are a helpful assistant. You have access to the following tools:\n\n{tools}\n\nPlease respond to the user's request. If a tool is needed, use it.",
        ),
        MessagesPlaceholder(variable_name="messages"),
    ]
)
