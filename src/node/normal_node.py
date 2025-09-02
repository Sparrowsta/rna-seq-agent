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
from langchain.tools import Tool

def create_normal_agent():
    """创建Normal节点的React Agent - 支持结构化输出"""
    # 使用共享的LLM实例
    llm = get_shared_llm()
    
    # 系统提示词 - 使用集中管理的prompt
    system_prompt = NORMAL_NODE_PROMPT
    
    tools = [
        # 核心项目信息中心工具
        Tool(
            name="get_project_overview",
            func=get_project_overview,
            description="项目全貌概览 - 一键查看项目完整状态和健康度。整合FASTQ数据、基因组状态、历史分析和资源评估。当用户询问'项目概览'、'项目状态'、'整体情况'时，立即调用此工具。"
        ),
        Tool(
            name="list_analysis_history",
            func=list_analysis_history,
            description="历史分析管理 - 浏览和管理已完成的分析。显示分析记录、结果大小、分析步骤和可复用配置。当用户询问'历史分析'、'分析记录'、'历史结果'时调用此工具。"
        ),
        
        # 详细信息查询工具（使用新的双模式工具）
        Tool(
            name="scan_fastq_files",
            func=lambda query="": scan_fastq_files(mode="normal", depth="detailed"),
            description="详细FASTQ文件分析 - 在整个项目目录递归扫描并列出所有可用的FASTQ文件。提供智能概览、统计信息、分析建议和详细样本信息。当用户询问'查看FASTQ文件'、'测序数据'、'数据文件'时调用此工具。"
        ),
        Tool(
            name="scan_genome_files", 
            func=lambda query="": scan_genome_files(mode="normal"),
            description="基因组配置查询 - 自动列出系统中所有可用的参考基因组。显示基因组版本、下载状态和文件大小。当用户询问'基因组'、'参考基因组'、'基因组信息'时调用此工具。"
        ),
        Tool(
            name="add_genome_config",
            func=add_genome_config,
            description="智能基因组配置 - 添加基因组配置。当用户说'添加基因组'并提供URL时调用此工具。直接传递用户的完整输入内容，工具内部会智能解析URL并提取基因组信息。"
        ),
        
        # 帮助工具
        Tool(
            name="get_help",
            func=get_help,
            description="系统功能帮助 - 显示Normal模式（项目信息中心）的完整功能列表和使用指南。当用户询问'帮助'、'功能'、'怎么用'时调用此工具。"
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
        
        # LangGraph的create_react_agent使用response_format时，结构化输出在result["structured_response"]中
        structured_response = result.get("structured_response")
        
        query_response = structured_response.query_response
        user_requirements = structured_response.user_requirements
        
        return {
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
            "messages": state.messages,
            "query_response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "status": "error"
        }