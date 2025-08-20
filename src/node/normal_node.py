import os
from typing import Dict, Any
from ..state import AgentState, NormalResponse
from ..tools import (
    query_fastq_files, 
    query_genome_info, 
    list_directory_tree,
    get_help
)
from ..core import get_shared_llm

from langgraph.prebuilt import create_react_agent
from langchain.tools import Tool

def create_normal_agent():
    """创建Normal节点的React Agent - 支持结构化输出"""
    # 使用共享的LLM实例
    llm = get_shared_llm()
    
    # 系统提示词 - 指导Agent行为和输出格式
    system_prompt = """你是RNA-seq智能分析助手的信息查询助手。你的任务是理解用户需求并调用合适的工具获取信息。

重要指导原则：
1. 根据用户的查询内容，选择最合适的工具
2. 调用工具后，将工具返回的完整结果作为你的最终回复
3. 不要对工具结果进行总结或解释，直接展示原始结果
4. 确保工具的输出信息完整传递给用户

可用工具说明：
- query_fastq_files: 当用户询问"FASTQ文件"、"测序数据"、"数据文件"时使用
- query_genome_info: 当用户询问"基因组"、"参考基因组"、"基因组信息"时使用  
- list_directory_tree: 当用户询问"目录结构"、"项目结构"、"文件结构"时使用
- search_ucsc_genomes: 当用户询问"搜索基因组"、"UCSC基因组"时使用
- get_help: 当用户询问"帮助"、"功能"、"使用方法"时使用

请直接调用工具并返回工具的完整输出结果。"""
    
    tools = [
        # 信息查询工具
        Tool(
            name="query_fastq_files",
            func=query_fastq_files,
            description="在整个项目目录递归扫描并列出所有可用的FASTQ文件。自动识别测序数据分布的各个目录，过滤掉已处理的文件。当用户询问'查看fastq文件'、'有什么数据文件'、'查看所有测序数据'时，立即调用此工具。"
        ),
        Tool(
            name="query_genome_info", 
            func=query_genome_info,
            description="自动列出系统中所有可用的参考基因组。无需用户提供具体信息，工具会自动显示支持的基因组版本和本地文件状态。当用户询问'基因组'、'参考基因组'、'支持什么基因组'时，立即调用此工具。"
        ),
        Tool(
            name="list_directory_tree",
            func=list_directory_tree,
            description="显示项目目录结构，帮助用户了解文件组织。当用户询问'目录结构'、'文件在哪里'、'查看项目结构'时调用此工具。"
        ),
        
        # 帮助工具
        Tool(
            name="get_help",
            func=get_help,
            description="显示系统功能帮助信息。当用户询问'帮助'、'功能'、'怎么用'时，立即调用此工具。"
        )
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
        
        structured_response = result.get("structured_response")
        
        if structured_response:
            normal_data = structured_response.dict()
            
            return {
                "messages": state.messages,
                "response": normal_data.get("query_response", ""),
                "status": "completed",
                "query_type": normal_data.get("query_type", "info"),
                "routing_decision": "normal",
                "query_response": normal_data.get("query_response", ""),
                "user_intent": normal_data.get("user_intent", ""),
                "suggested_actions": normal_data.get("suggested_actions", [])
            }
        else:
            raise Exception("Agent未返回预期的结构化响应")
        
    except Exception as e:
        return {
            "messages": state.messages,
            "query_response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "status": "error",
            "routing_decision": "normal"
        }