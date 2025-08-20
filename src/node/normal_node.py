import os
from typing import Dict, Any
from ..state import AgentState
from ..tools import query_fastq_files, query_genome_info, get_help
from ..core import get_shared_llm

from langgraph.prebuilt import create_react_agent
from langchain.tools import Tool

def create_normal_agent():
    """创建Normal节点的React Agent - 支持结构化输出"""
    # 使用共享的LLM实例
    llm = get_shared_llm()
    
    tools = [
        Tool(
            name="query_fastq_files",
            func=query_fastq_files,
            description="自动扫描并列出系统中所有可用的FASTQ文件。无需用户提供路径，工具会自动发现data/fastq目录下的所有文件。当用户询问'查看fastq文件'、'有什么数据文件'等时，立即调用此工具。"
        ),
        Tool(
            name="query_genome_info", 
            func=query_genome_info,
            description="自动列出系统中所有可用的参考基因组。无需用户提供具体信息，工具会自动显示支持的基因组版本。当用户询问'基因组'、'参考基因组'、'支持什么基因组'时，立即调用此工具。"
        ),
        Tool(
            name="get_help",
            func=get_help,
            description="显示系统功能帮助信息。当用户询问'帮助'、'功能'、'怎么用'时，立即调用此工具。"
        )
    ]
    
    # 使用LangGraph预构件，支持Pydantic结构化输出
    agent = create_react_agent(
        model=llm,
        tools=tools,
        response_format=AgentState  # 关键：使用Pydantic模型作为响应格式
    )
    return agent

async def normal_node(state: AgentState) -> Dict[str, Any]:
    """Normal节点 - 使用LangGraph React Agent预构件处理用户查询"""
    
    try:
        agent_executor = create_normal_agent()
        result = await agent_executor.ainvoke({
            "messages": state.messages
        })
        
        # LangGraph会自动生成符合AgentState格式的结构化响应
        structured_response = result.get("structured_response")
        
        if structured_response:
            # 返回结构化的Pydantic模型数据
            return structured_response.dict()
        else:
            # 后备方案：手动构造响应
            output = result.get("output", "处理完成")
            
            return {
                "query_type": "tool_executed",
                "routing_decision": "normal", 
                "query_response": output,
                "user_intent": "React Agent处理完成",
                "suggested_actions": ["继续查询", "开始分析"],
                "response": output,
                "status": "completed"
            }
        
    except Exception as e:
        return {
            "query_type": "error",
            "routing_decision": "normal",
            "query_response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "user_intent": "系统错误", 
            "suggested_actions": ["重试"],
            "response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "status": "error"
        }