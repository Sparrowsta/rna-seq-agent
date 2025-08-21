import os
from typing import Dict, Any
from ..state import AgentState, NormalResponse
from ..tools import (
    query_fastq_files, 
    query_genome_info, 
    get_help,
    add_genome_config,
    get_project_overview,
    list_analysis_history
)
from ..core import get_shared_llm

from langgraph.prebuilt import create_react_agent
from langchain.tools import Tool

def create_normal_agent():
    """创建Normal节点的React Agent - 支持结构化输出"""
    # 使用共享的LLM实例
    llm = get_shared_llm()
    
    # 系统提示词 - 指导Agent行为和输出格式
    system_prompt = """你是RNA-seq智能分析助手的项目信息中心。你的任务是理解用户需求并调用合适的工具获取信息。

重要指导原则：
1. 根据用户的查询内容，选择最合适的工具并主动调用
2. 调用工具后，将工具返回的完整结果作为你的最终回复
3. 不要对工具结果进行总结或解释，直接展示原始结果
4. 确保工具的输出信息完整传递给用户

核心项目工具：
- get_project_overview: 当用户询问"项目概览"、"项目状态"、"整体情况"时使用
- list_analysis_history: 当用户询问"历史分析"、"分析记录"、"历史结果"时使用

详细查询工具：
- query_fastq_files: 当用户询问"FASTQ文件"、"测序数据"、"数据文件"时使用
- query_genome_info: 当用户询问"基因组"、"参考基因组"、"基因组信息"时使用  
- add_genome_config: 当用户说"添加基因组"并提供URL时，直接传递完整的用户输入
- get_help: 当用户询问"帮助"、"功能"、"使用方法"时使用

请直接调用工具并返回工具的完整输出结果。"""
    
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
        
        # 详细信息查询工具
        Tool(
            name="query_fastq_files",
            func=query_fastq_files,
            description="详细FASTQ文件分析 - 在整个项目目录递归扫描并列出所有可用的FASTQ文件。提供智能概览、统计信息、分析建议和详细样本信息。当用户询问'查看FASTQ文件'、'测序数据'、'数据文件'时调用此工具。"
        ),
        Tool(
            name="query_genome_info", 
            func=query_genome_info,
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
        print(f"🔍 Normal节点开始处理，最新消息: {state.messages[-1] if state.messages else '无消息'}")
        
        agent_executor = create_normal_agent()
        messages_input = {"messages": state.messages}
        
        print(f"📨 传入Agent的消息: {[getattr(msg, 'content', str(msg)) for msg in state.messages]}")
        
        result = await agent_executor.ainvoke(messages_input)
        
        print(f"📋 Agent返回结果: {result}")
        
        structured_response = result.get("structured_response")
        
        if structured_response:
            normal_data = structured_response.dict()
            print(f"✅ 结构化响应: {normal_data}")
            
            return {
                "messages": state.messages,
                "query_response": normal_data.get("query_response", ""),
                "status": "completed",
                "query_type": normal_data.get("query_type", "info"),
                "user_intent": normal_data.get("user_intent", ""),
                "suggested_actions": normal_data.get("suggested_actions", [])
            }
        else:
            print("❌ Agent未返回预期的结构化响应")
            raise Exception("Agent未返回预期的结构化响应")
        
    except Exception as e:
        print(f"❌ Normal节点处理出错: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return {
            "messages": state.messages,
            "query_response": f"抱歉，处理您的请求时出现错误: {str(e)}",
            "status": "error"
        }