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
    system_prompt = """你是RNA-seq智能分析助手的项目信息中心。你的核心任务是：
1. 理解用户需求并调用合适的工具获取信息
2. **重要：智能识别用户的配置需求并输出结构化配置**

🔍 **双重输出策略**：
- query_response: 工具调用的完整结果
- user_requirements: 用户配置需求的结构化字典格式（类似nextflow_config）

📝 **结构化配置需求提取规则**：
根据用户输入识别并输出标准化的配置字典：

**基因组选择识别**：
- "使用hg19/选择hg19/hg19基因组" → {"genome_version": "hg19", "species": "human"}
- "使用hg38/选择hg38/hg38基因组" → {"genome_version": "hg38", "species": "human"}  
- "使用mm10/选择mm10/mm10基因组" → {"genome_version": "mm10", "species": "mouse"}
- "使用mm39/选择mm39/mm39基因组" → {"genome_version": "mm39", "species": "mouse"}

**工具选择识别**：
- "用fastp/fastp质控/选择fastp" → {"qc_tool": "fastp"}
- "用cutadapt/cutadapt质控" → {"qc_tool": "cutadapt"}
- "用STAR/STAR比对/选择STAR" → {"align_tool": "star"}
- "用hisat2/hisat2比对" → {"align_tool": "hisat2"}  
- "用featureCounts/featureCounts定量" → {"quant_tool": "featurecounts"}
- "用htseq/htseq定量" → {"quant_tool": "htseq"}

**分析类型识别**：
- "差异表达/差异基因/找差异基因" → {"analysis_type": "differential_expression"}
- "质量控制/质控分析/数据质控" → {"analysis_type": "quality_control"}

**其他配置识别**：
- "双端测序/paired-end/PE数据" → {"paired_end": true}
- "单端测序/single-end/SE数据" → {"paired_end": false}

💡 **处理示例**：

用户输入："使用hg19进行差异基因分析"
→ query_response: "[调用相关工具的结果]"
→ user_requirements: {"genome_version": "hg19", "species": "human", "analysis_type": "differential_expression"}

用户输入："用STAR和fastp分析RNA数据"
→ query_response: "[调用相关工具的结果]"  
→ user_requirements: {"qc_tool": "fastp", "align_tool": "star"}

用户输入："查看基因组信息"（只是查询，无配置意图）
→ query_response: "[基因组查询工具的结果]"
→ user_requirements: {}

核心项目工具：
- get_project_overview: 当用户询问"项目概览"、"项目状态"、"整体情况"时使用
- list_analysis_history: 当用户询问"历史分析"、"分析记录"、"历史结果"时使用

详细查询工具：
- query_fastq_files: 当用户询问"FASTQ文件"、"测序数据"、"数据文件"时使用
- query_genome_info: 当用户询问"基因组"、"参考基因组"、"基因组信息"时使用  
- add_genome_config: 当用户说"添加基因组"并提供URL时，直接传递完整的用户输入
- get_help: 当用户询问"帮助"、"功能"、"使用方法"时使用

请调用工具并返回完整的结构化输出，包括工具结果和结构化的配置需求。"""
    
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