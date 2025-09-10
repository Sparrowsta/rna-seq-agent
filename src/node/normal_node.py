from typing import Dict, Any
from ..state import AgentState, NormalResponse
from ..tools import (
    scan_fastq_files,
    scan_genome_files, 
    get_help,
    add_genome_config,
    get_project_overview,
    list_analysis_history
)
from ..core import get_shared_llm
from ..prompts import NORMAL_NODE_PROMPT

from langgraph.prebuilt import create_react_agent

def create_normal_agent():
    """创建Normal节点的React Agent - 支持结构化输出"""
    # 使用共享的LLM实例
    llm = get_shared_llm()
    
    # 系统提示词 - 使用集中管理的prompt
    system_prompt = NORMAL_NODE_PROMPT
    
    # 直接使用@tool装饰的函数，无需Tool包装
    tools = [
        get_project_overview,
        list_analysis_history,
        scan_fastq_files,
        scan_genome_files,
        add_genome_config,
        get_help
    ]
    
    # 使用LangGraph预构件，使用精简的响应格式
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=system_prompt,  # 添加系统提示词
        response_format=NormalResponse  # 使用精简的响应格式
    )
    return agent

async def normal_node(state: AgentState) -> Dict[str, Any]:
    """Normal节点 - 使用LangGraph React Agent预构件处理用户查询"""
    
    try:
        agent_executor = create_normal_agent()
        messages_input = {"messages": state.messages}
        
        result = await agent_executor.ainvoke(messages_input)
        
        # LangGraph的create_react_agent使用response_format时，结构化输出在result["structured_response"]中
        structured_response = result.get("structured_response")
        
        query_response = structured_response.query_response
        user_requirements = structured_response.user_requirements
        
        return {
            "success": True,
            "messages": result.get("messages", state.messages),
            "query_response": query_response,
            "user_requirements": user_requirements,
            "status": "normal"
        }
        
    except Exception as e:
        print(f"❌ Normal节点处理出错: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return {
            "success": False,
            "messages": state.messages,
            "query_response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "status": "error"
        }
