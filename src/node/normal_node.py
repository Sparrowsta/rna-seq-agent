from typing import Dict, Any
from ..state import NormalNodeState

async def normal_node(state: NormalNodeState) -> Dict[str, Any]:
    """Normal节点 - 用户交互和信息查询"""
    print(f"💬 处理用户请求...")
    print(f"   用户输入: {state['input']}")
    
    # TODO: 实现用户输入分析和分类逻辑
    # TODO: 实现信息查询功能（FASTQ文件、基因组列表等）
    # TODO: 实现路由决策逻辑（normal/plan）
    # TODO: 实现帮助信息和引导功能
    
    return {
        "query_type": "info",  # 临时默认值
        "routing_decision": "normal",  # 临时默认保持在normal
        "query_response": "用户交互功能待实现",
        "response": "已接收用户输入，处理功能开发中",
        "status": "normal"
    }