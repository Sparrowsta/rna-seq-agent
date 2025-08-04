import os
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from .tools import list_directory, read_genomes_config, update_genomes_config, ReadGenomeConfigArgs, UpdateGenomeConfigArgs

tools = [list_directory, read_genomes_config, update_genomes_config]
llm = ChatOpenAI(
    model=os.environ.get("OPENAI_MODEL_NAME"),
    api_key=os.environ.get("OPENAI_API_KEY"),
    base_url=os.environ.get("OPENAI_API_BASE")
)
llm_with_tools = llm.bind_tools(tools)

prompt = ChatPromptTemplate.from_messages(
    [
        (
            "system",
            """你是一位专家级的生物信息学家助手。你的主要职责是帮助用户管理基因组配置文件。

**核心工作流程指南:**

1.  **`read_genomes_config`**: 使用此工具可以读取当前 `genomes.json` 配置文件的内容。这是了解现有配置的唯一方法。
2.  **`update_genomes_config`**: 使用此工具可以更新 `genomes.json` 配置文件。它会完全覆盖现有文件，因此请务必谨慎使用。你的输入是一个json格式内容

请根据用户的请求，并严格遵循上述指南来执行任务。""",
        ),
        MessagesPlaceholder(variable_name="messages"),
    ]
)
