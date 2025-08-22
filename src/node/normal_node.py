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
2. **重要：同时智能识别用户的配置意图，在config_updates字段中返回相应的Nextflow参数**

🔍 **双重处理策略**：
- 如果用户询问信息（如"查看基因组"），调用相关工具
- 如果用户表达配置意图（如"使用hg19"），在config_updates中设置相应参数
- 两者可以同时进行

⚙️ **配置识别智能规则**：

**基因组选择识别**：
- 当用户提到"使用/用/选择/基于 + 基因组名称"时，设置genome_version
- 常见基因组：hg38/hg19(human), mm39/mm10/mm9(mouse), danRer11(zebrafish), xenLae2(xenopus), ce11(worm)
- 示例：
  * "使用hg19" → {"genome_version": "hg19", "species": "human"}
  * "用mm10" → {"genome_version": "mm10", "species": "mouse"}
  * "选择hg38基因组" → {"genome_version": "hg38", "species": "human"}

**工具选择识别**：
- QC工具：fastp, cutadapt → {"qc_tool": "工具名"}
- 比对工具：star, hisat2 → {"align_tool": "工具名"}  
- 定量工具：featurecounts, htseq → {"quant_tool": "工具名"}

**分析类型识别**：
- "差异表达/差异基因" → {"analysis_type": "differential_expression"}
- "质量控制/质控" → {"analysis_type": "quality_control"}

💡 **处理示例**：
用户输入："使用hg19"
- 调用query_genome_info工具查看基因组状态
- 同时在config_updates中设置：{"genome_version": "hg19", "species": "human"}

核心项目工具：
- get_project_overview: 当用户询问"项目概览"、"项目状态"、"整体情况"时使用
- list_analysis_history: 当用户询问"历史分析"、"分析记录"、"历史结果"时使用

详细查询工具：
- query_fastq_files: 当用户询问"FASTQ文件"、"测序数据"、"数据文件"时使用
- query_genome_info: 当用户询问"基因组"、"参考基因组"、"基因组信息"时使用  
- add_genome_config: 当用户说"添加基因组"并提供URL时，直接传递完整的用户输入
- get_help: 当用户询问"帮助"、"功能"、"使用方法"时使用

请直接调用工具并返回工具的完整输出结果，同时识别配置意图并更新config_updates字段。"""
    
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
        
        print(f"📋 Agent返回结果类型: {type(result)}")
        print(f"📋 Agent返回结果: {result}")
        
        # 探索返回结构
        if isinstance(result, dict):
            print(f"🔍 字典keys: {list(result.keys())}")
            for key, value in result.items():
                print(f"   {key}: {type(value)} - {str(value)[:100]}...")
        
        structured_response = result.get("structured_response")
        print(f"🎯 structured_response: {structured_response}")
        print(f"🎯 structured_response类型: {type(structured_response)}")
        
        if structured_response:
            print(f"✅ 结构化响应: {structured_response}")
            
            # 直接使用Pydantic模型属性
            config_updates = structured_response.config_updates or {}
            query_response = structured_response.query_response or ""
            
            # nextflow_config已在state中初始化，无需检查None
            updated_nextflow_config = state.nextflow_config.copy()
            
            if config_updates:
                print(f"🔧 检测到配置更新: {config_updates}")
                updated_nextflow_config.update(config_updates)
                print(f"📝 更新后的nextflow_config: {updated_nextflow_config}")
            
            return {
                "messages": state.messages,
                "query_response": query_response,
                "status": "normal",
                "nextflow_config": updated_nextflow_config
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