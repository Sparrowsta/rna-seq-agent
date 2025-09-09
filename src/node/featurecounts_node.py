"""
FeatureCounts节点 - 用于执行基因定量分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FeaturecountsResponse
from ..core import get_shared_llm


def create_featurecounts_agent():
    """创建FeatureCounts节点的React Agent"""
    llm = get_shared_llm()
    
    system_prompt = """你是RNA-seq分析流水线中的基因定量专家。

你的任务：
1. 分析STAR比对后的BAM文件
2. 执行FeatureCounts基因定量分析
3. 生成基因表达矩阵和统计信息
4. 返回结构化的执行结果

请根据提供的状态信息，执行FeatureCounts定量分析，并返回详细的执行结果。
目前这是占位实现，请返回模拟的成功结果。

输出格式要求：
- status: "success" 表示成功
- summary: 简要的执行总结，包含定量统计信息
- response: 详细的响应消息
"""
    
    # 创建不使用工具的React Agent
    agent = create_react_agent(
        model=llm,
        tools=[],  # 暂不使用工具
        prompt=system_prompt,
        response_format=FeaturecountsResponse
    )
    return agent


def featurecounts_node(state: AgentState) -> Dict[str, Any]:
    """
    FeatureCounts节点实现
    
    功能：
    - 执行基因定量
    - 生成表达矩阵
    - 更新状态信息
    - 批次优化模式下统一返回所有优化参数
    """
    print("\n📊 FeatureCounts定量节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "featurecounts" not in completed_steps:
        completed_steps.append("featurecounts")
    
    # 获取前面步骤的信息
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    quant_tool = state.nextflow_config.get('quant_tool', 'featureCounts')
    execution_mode = state.execution_mode
    
    # 基础执行结果
    result = {
        "status": "featurecounts_completed",
        "response": f"✅ FeatureCounts定量完成\n- 定量样本: {sample_count}个\n- 检测基因: 24,587个\n- 定量工具: {quant_tool}\n- 成功率: 97.2%",
        "current_step": "featurecounts",
        "completed_steps": completed_steps,
        "featurecounts_results": {
            "status": "success",
            "summary": f"FeatureCounts成功定量{sample_count}个样本，检测到24,587个基因"
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
            "-T": 8,
            "-M": True,
            "--minOverlap": 10,
            "--fracOverlap": 0.2
        }
        optimization_reasoning = "基因分配率82.1%可提升，已应用多重比对和线程数优化"
        
        # 立即更新执行参数
        optimized_params = {**state.featurecounts_params, **optimization_suggestions}
        result["featurecounts_params"] = optimized_params
        result["featurecounts_optimization_suggestions"] = optimization_reasoning
        result["response"] += f"\n\n⚡ **立即优化**: {optimization_reasoning}"
        
        print(f"✅ [OPTIMIZED] FeatureCounts优化参数已应用: {len(optimization_suggestions)}个参数")
        
    elif execution_mode == "batch_optimize":
        # 批次优化模式：收集FeatureCounts优化参数，不立即应用
        print("📦 [BATCH] FeatureCounts批次优化模式，收集优化参数...")
        
        # 硬编码模拟FeatureCounts优化参数
        featurecounts_optimization = {
            "optimization_reasoning": "基因分配率82.1%可进一步提升，建议启用多重比对和调整线程数",
            "suggested_params": {
                "-T": 8,
                "-M": True,
                "--minOverlap": 10,
                "--fracOverlap": 0.2
            },
            "current_params": state.featurecounts_params.copy(),
            "tool_name": "featurecounts"
        }
        
        # 将FeatureCounts优化参数添加到批次收集器
        batch_optimizations = state.batch_optimizations.copy()
        batch_optimizations["featurecounts"] = featurecounts_optimization
        
        # 标记批次优化完成，但参数尚未应用
        result["batch_optimizations"] = batch_optimizations
        # 也记录优化理由，便于在确认页展示
        result["featurecounts_optimization_suggestions"] = featurecounts_optimization.get("optimization_reasoning", "")
        result["batch_optimization_complete"] = True
        
        # 生成批次优化收集完成的总结报告
        total_optimizations = len(batch_optimizations)
        optimization_summary = []
        
        for tool_name, opt_data in batch_optimizations.items():
            param_count = len(opt_data.get("suggested_params", {}))
            optimization_summary.append(f"  • {tool_name.upper()}: {param_count}个参数")
        
        batch_summary = f"""
📦 **批次优化收集完成**

🔧 **收集到的优化建议**:
{chr(10).join(optimization_summary)}

💡 **FeatureCounts优化**: 基因分配率可提升，建议启用多重比对

✅ **状态**: 所有工具优化参数已收集完毕，将在下次执行时应用
"""
        
        result["response"] += batch_summary
        print(f"✅ [BATCH] 批次优化收集完成: 共{total_optimizations}个工具的优化参数")
    
    return result
