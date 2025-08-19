from typing import Dict, Any
from langgraph.graph import StateGraph, START, END
from .state import NormalNodeState
from .node.normal_node import normal_node
from .node.plan_node import plan_node
from .node.detect_node import detect_node
from .node.prepare_node import prepare_node
from .node.replan_node import replan_node
from .node.user_confirm_node import user_confirm_node
from .node.execute_node import execute_node
from .route import route_from_normal, should_continue, route_after_confirm

def create_agent():
    """创建LangGraph Agent - 完整的Plan-and-Execute架构"""
    
    # 创建状态图
    workflow = StateGraph(NormalNodeState)
    
    # 添加七个专业化节点
    workflow.add_node("normal", normal_node)
    workflow.add_node("plan", plan_node)
    workflow.add_node("detect", detect_node)
    workflow.add_node("prepare", prepare_node)
    workflow.add_node("user_confirm", user_confirm_node)
    workflow.add_node("replan", replan_node)
    workflow.add_node("execute", execute_node)
    
    # 入口点：从START进入Normal节点
    workflow.add_edge(START, "normal")
    
    # Normal节点的条件路由
    workflow.add_conditional_edges(
        "normal",
        route_from_normal,
        {
            "normal": "normal",      # 继续交互
            "plan": "plan"           # 进入分析流程
        }
    )
    
    # 主流程: Normal → Plan → Detect → Prepare → UserConfirm
    workflow.add_edge("plan", "detect")
    workflow.add_edge("detect", "prepare")
    workflow.add_edge("prepare", "user_confirm")
    
    # 用户确认后的条件路由
    workflow.add_conditional_edges(
        "user_confirm",
        route_after_confirm,
        {
            "execute": "execute",    # 执行分析 - 进入execute节点
            "modify": "replan",      # 修改配置 - 进入replan流程
            "cancel": END            # 取消分析 - 结束图流程
        }
    )
    
    # Execute节点执行完成后直接结束
    workflow.add_edge("execute", END)
    
    # Replan流程: Replan → Detect → Prepare → UserConfirm
    workflow.add_edge("replan", "detect")
    
    # 编译图
    app = workflow.compile()
    
    print("🤖 完整Plan-and-Execute Agent已创建 (Normal → Plan → Detect → Prepare → UserConfirm + Replan)")
    return app