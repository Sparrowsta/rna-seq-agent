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
    execution_mode = state.execution_mode
    
    # 基础执行结果
    result = {
        "status": "star_completed",
        "response": f"✅ STAR比对完成\n- 比对样本: {sample_count}个\n- 物种: {species}\n- 比对率: 88.5%\n- 唯一比对: 82.3%",
        "current_step": "star",
        "completed_steps": completed_steps,
        "star_results": {
            "status": "success",
            "summary": f"STAR成功比对{sample_count}个样本，平均比对率88.5%"
        }
    }
    
    # 根据执行模式处理优化逻辑
    if execution_mode == "single":
        # 单次执行模式：不生成任何优化参数
        print("🚀 [SINGLE] 单次执行模式，直接完成")
        result["response"] += "\n\n🚀 **单次执行**: 任务完成，无优化处理"
        
    elif execution_mode == "optimized":
        # 精细优化模式：立即应用优化参数
        print("⚡ [OPTIMIZED] 精细优化模式，立即应用优化...")
        
        # 硬编码模拟优化参数
        optimization_suggestions = {
            "--runThreadN": 16,
            "--outFilterMultimapNmax": 20,
            "--outFilterMismatchNmax": 2,
            "--alignIntronMax": 1000000
        }
        optimization_reasoning = "比对率88.5%略低，已应用多重比对和线程数优化"
        
        # 立即更新执行参数
        optimized_params = {**state.star_params, **optimization_suggestions}
        result["star_params"] = optimized_params
        result["star_optimization_suggestions"] = optimization_reasoning
        result["response"] += f"\n\n⚡ **立即优化**: {optimization_reasoning}"
        
        print(f"✅ [OPTIMIZED] STAR优化参数已应用: {len(optimization_suggestions)}个参数")
        
    elif execution_mode == "batch_optimize":
        # 批次优化模式：收集优化参数
        print("📦 [BATCH] STAR批次优化模式，收集优化参数...")
        
        # 硬编码模拟STAR优化参数
        star_optimization = {
            "optimization_reasoning": "比对率88.5%略低于最佳水平，建议调整多重比对参数和线程数以提升性能",
            "suggested_params": {
                "--runThreadN": 16,
                "--outFilterMultimapNmax": 20,
                "--outFilterMismatchNmax": 2,
                "--alignIntronMax": 1000000
            },
            "current_params": state.star_params.copy(),
            "tool_name": "star"
        }
        
        # 将优化参数添加到批次收集器（继承之前的收集结果）
        batch_optimizations = state.batch_optimizations.copy()
        batch_optimizations["star"] = star_optimization
        
        result["batch_optimizations"] = batch_optimizations
        # 也记录优化理由，便于在确认页展示
        result["star_optimization_suggestions"] = star_optimization.get("optimization_reasoning", "")
        result["response"] += "\n\n📦 **STAR优化参数已收集**: 比对率偏低，建议调整多重比对参数"
        
        print(f"✅ [BATCH] STAR优化参数收集完成: {len(star_optimization['suggested_params'])}个参数")
    
    return result
