from typing import Dict, Any, List
from ..state import AgentState, DetectResponse
from langgraph.prebuilt import create_react_agent
from langchain.tools import Tool
from ..tools import (
    analyze_fastq_data,
    assess_system_readiness,
    verify_genome_setup,
    check_fastp_availability,
    check_star_availability,
    check_featurecounts_availability
)
from ..core import get_shared_llm


def create_detection_agent():
    """创建检测执行Agent"""
    llm = get_shared_llm()
    
    # 直接定义检测工具列表
    tools = [
        Tool(
            name="analyze_fastq_data",
            func=analyze_fastq_data,
            description="扫描和分析FASTQ文件。收集项目中所有FASTQ文件的信息，包括文件大小、样本配对关系、测序类型等。"
        ),
        Tool(
            name="assess_system_readiness", 
            func=assess_system_readiness,
            description="检测系统硬件资源。评估CPU核心数、内存容量、磁盘空间、系统负载等硬件信息。"
        ),
        Tool(
            name="verify_genome_setup",
            func=verify_genome_setup, 
            description="验证基因组文件配置。检查已配置基因组的FASTA文件、GTF文件、STAR索引文件的存在性和完整性。"
        ),
        Tool(
            name="check_fastp_availability",
            func=check_fastp_availability,
            description="检测fastp质控工具的可用性。在micromamba环境中测试fastp命令是否可执行。"
        ),
        Tool(
            name="check_star_availability", 
            func=check_star_availability,
            description="检测STAR比对工具的可用性。在micromamba环境中测试STAR命令是否可执行。"
        ),
        Tool(
            name="check_featurecounts_availability",
            func=check_featurecounts_availability,
            description="检测featureCounts定量工具的可用性。在micromamba环境中测试featureCounts工具是否可执行。"
        )
    ]
    
    system_prompt = """你是检测执行专家。你的任务是根据计划列表，智能执行检测任务并收集结果。

执行原则：
1. 按计划列表中的任务名称，依次调用对应的检测工具
2. 对于每个任务，只调用一次相应的工具
3. 如果检测失败，记录错误但继续执行其他任务
4. 收集所有检测结果，整合成统一的数据结构

返回格式示例：
{
  "query_results": {"检测结果按任务整理": "工具已优化输出格式"},
  "query_summary": "检测完成：基因组hg19可用，工具就绪，发现6个FASTQ文件"
}

可用的检测工具：
- analyze_fastq_data: 分析FASTQ文件
- assess_system_readiness: 检测系统资源
- verify_genome_setup: 验证基因组配置
- check_fastp_availability: 检测fastp工具
- check_star_availability: 检测STAR工具
- check_featurecounts_availability: 检测featureCounts工具

请按照计划列表执行检测，并返回JSON格式的结果。"""
    
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=system_prompt,
        response_format=DetectResponse
    )
    return agent


async def detect_node(state: AgentState) -> Dict[str, Any]:
    """Detect节点 - 使用LangGraph React Agent执行检测任务"""
    print(f"🔍 Detect节点开始处理，计划任务: {state.plan}")
    
    if not state.plan:
        return {
            "query_summary": "没有检测任务需要执行",
            "status": "prepare",
            "query_results": {}
        }
    
    try:
        # 构建执行计划的消息
        task_list = "\n".join([f"{i+1}. {task}" for i, task in enumerate(state.plan)])
        messages = [{
            "role": "user", 
            "content": f"请按顺序执行以下检测任务：\n\n{task_list}\n\n请调用相应的工具并收集所有检测数据。"
        }]
        
        agent_executor = create_detection_agent()
        messages_input = {"messages": messages}
        
        print(f"🤖 开始执行检测任务")
        result = await agent_executor.ainvoke(messages_input)
        
        print(f"📋 Agent返回结果类型: {type(result)}")
        print(f"📋 Agent返回结果: {result}")
        
        structured_response = result.get("structured_response")
        print(f"🎯 structured_response: {structured_response}")
        
        if structured_response:
            print(f"✅ 检测成功，准备返回数据")
            
            query_results = structured_response.query_results or {}
            query_summary = structured_response.query_summary or "检测完成"
            
            print(f"📊 query_results keys: {list(query_results.keys())}")
            print(f"📝 query_summary: {query_summary}")
            
            return {
                "query_summary": query_summary,
                "status": "prepare", 
                "query_results": query_results
            }
        else:
            print("❌ Agent未返回预期的结构化响应")
            raise Exception("Agent未返回预期的结构化响应")
            
    except Exception as e:
        print(f"❌ 检测执行失败: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return {
            "query_summary": f"检测失败: {str(e)}",
            "status": "error",
            "query_results": {}
        }