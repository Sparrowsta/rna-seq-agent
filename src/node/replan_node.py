from typing import Dict, Any
from ..state import ReplanNodeState

async def replan_node(state: ReplanNodeState) -> Dict[str, Any]:
    """重新计划节点 - 分析用户修改请求并决定路由"""
    print(f"🔄 重新计划中...")
    print(f"   用户修改请求: {state['input']}")
    
    # TODO: 实现LLM驱动的修改意图分析
    # TODO: 实现智能路由决策逻辑
    # 参考docs/RNA-seq_Agent_开发计划.md中的ReplanResponse模型
    
    # 临时基础实现 - 统一路由到detect节点
    return {
        "user_modification_input": state.get('input', ''),
        "modification_intent": {},
        "modification_mode": "incremental",
        "routing_decision": "detect",
        "routing_reason": "用户修改需要重新检测系统状态",
        "response": "已接收修改请求，将重新检测系统状态",
        "status": "replanning"
    }