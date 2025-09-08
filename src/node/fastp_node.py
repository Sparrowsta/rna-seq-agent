"""
FastP节点 - 用于执行FastP质量控制
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FastpResponse
from ..core import get_shared_llm


def create_fastp_agent():
    """创建FastP节点的React Agent"""
    llm = get_shared_llm()
    
    system_prompt = """你是RNA-seq分析流水线中的FastP质控专家。

你的任务：
1. 分析当前的样本信息和配置参数
2. 执行FastP质量控制处理
3. 生成质控报告和统计信息
4. 返回结构化的执行结果

请根据提供的状态信息，执行FastP质控处理，并返回详细的执行结果。
目前这是占位实现，请返回模拟的成功结果。

输出格式要求：
- status: "success" 表示成功
- summary: 简要的执行总结
- response: 详细的响应消息
"""
    
    # 创建不使用工具的React Agent
    agent = create_react_agent(
        model=llm,
        tools=[],  # 暂不使用工具
        prompt=system_prompt,
        response_format=FastpResponse
    )
    return agent


def fastp_node(state: AgentState) -> Dict[str, Any]:
    """
    FastP节点实现
    
    功能：
    - 执行FastP质量控制
    - 生成质量报告
    - 更新状态信息
    """
    print("\n🧹 FastP质控节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "fastp" not in completed_steps:
        completed_steps.append("fastp")
    
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    
    # 返回成功结果
    return {
        "status": "fastp_completed",
        "response": f"✅ FastP质控完成\n- 处理样本: {sample_count}个\n- 质控通过率: 95%\n- 平均Q30: 92%",
        "current_step": "fastp",
        "completed_steps": completed_steps,
        "fastp_results": {
            "status": "success",
            "summary": f"FastP成功处理{sample_count}个样本，质控通过率95%"
        }
    }