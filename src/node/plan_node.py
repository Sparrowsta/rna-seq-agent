from typing import Dict, Any
from ..state import PlanNodeState

async def plan_node(state: PlanNodeState) -> Dict[str, Any]:
    """计划制定节点"""
    print(f"🎯 制定计划中...")
    print(f"   用户输入: {state['input']}")
    
    # TODO: 实现计划制定逻辑
    # 这里需要基于用户输入生成分析步骤计划
    # 参考docs/plan-and-execute.ipynb中的planner实现
    
    return {
        "plan": ["检测FASTQ文件", "配置基因组参数", "准备执行流水线"],
        "analysis_intent": "RNA-seq差异基因表达分析",
        "response": "分析计划已制定完成",
        "status": "planning"
    }
