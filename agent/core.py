import os
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

llm = ChatGoogleGenerativeAI(
    model=os.environ.get("GEMINI_MODEL_NAME"),
    google_api_key=os.environ.get("GOOGLE_API_KEY")
)

prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            "You are a helpful assistant.",
        ),
        MessagesPlaceholder(variable_name="messages"),
    ]
)
