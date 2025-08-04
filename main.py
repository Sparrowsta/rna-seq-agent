from dotenv import load_dotenv
load_dotenv(dotenv_path='config/.env')

import os
from agent.graph import agent_executor
from langchain_core.messages import SystemMessage


if __name__ == "__main__":
    agent_executor.invoke({})