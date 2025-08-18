from typing import Dict, Any
from ..state import PlanExecuteState

async def execute_node(state: PlanExecuteState) -> Dict[str, Any]:
    """执行节点"""
    plan = state.get("plan", [])
    if not plan:
        return {"response": "无计划可执行"}
    
    current_step = plan[0]
    print(f"⚡ 执行步骤: {current_step}")
    
    # TODO(human): 实现步骤执行逻辑
    # 这里需要实际调用RNA-seq分析工具
    # 集成现有的Nextflow流水线
    
    return {
        "past_steps": [(current_step, f"步骤 '{current_step}' 执行完成")],
        "plan": plan[1:]  # 移除已完成的步骤
    }