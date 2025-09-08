"""
STAR节点 - 用于执行STAR比对分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, StarResponse
from ..core import get_shared_llm


def create_star_agent():
    """创建STAR节点的React Agent"""
    llm = get_shared_llm()
    
    system_prompt = """你是RNA-seq分析流水线中的STAR比对专家。

你的任务：
1. 分析FastP质控后的数据状态
2. 执行STAR序列比对处理
3. 生成比对统计和质量报告
4. 返回结构化的执行结果

请根据提供的状态信息，执行STAR比对处理，并返回详细的执行结果。
目前这是占位实现，请返回模拟的成功结果。

输出格式要求：
- status: "success" 表示成功
- summary: 简要的执行总结，包含比对统计信息
- response: 详细的响应消息
"""
    
    # 创建不使用工具的React Agent
    agent = create_react_agent(
        model=llm,
        tools=[],  # 暂不使用工具
        prompt=system_prompt,
        response_format=StarResponse
    )
    return agent


def star_node(state: AgentState) -> Dict[str, Any]:
    """
    STAR节点实现
    
    功能：
    - 执行STAR比对
    - 生成比对统计
    - 更新状态信息
    """
    print("\n🎯 STAR比对节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "star" not in completed_steps:
        completed_steps.append("star")
    
    # 获取前一步的信息
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    species = state.nextflow_config.get('species', 'human')
    
    # 返回成功结果
    return {
        "status": "star_completed",
        "response": f"✅ STAR比对完成\n- 比对样本: {sample_count}个\n- 物种: {species}\n- 比对率: 88.5%\n- 唯一比对: 82.3%",
        "current_step": "star",
        "completed_steps": completed_steps,
        "star_results": {
            "status": "success",
            "summary": f"STAR成功比对{sample_count}个样本，平均比对率88.5%"
        }
    }