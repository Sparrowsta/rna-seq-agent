from typing import Dict, Any
from langgraph.graph import StateGraph, START, END
from .state import PlanExecuteState
from .node.plan_node import plan_node
from .node.execute_node import execute_node
from .route import should_continue

def create_agent():
    """创建LangGraph Agent"""
    
    # 创建状态图
    workflow = StateGraph(PlanExecuteState)
    
    # 添加节点
    workflow.add_node("plan", plan_node)
    workflow.add_node("execute", execute_node)
    
    # 添加边
    workflow.add_edge(START, "plan")
    workflow.add_edge("plan", "execute")
    workflow.add_conditional_edges("execute", should_continue, ["execute", END])
    
    # 编译图
    app = workflow.compile()
    
    print("🤖 LangGraph Agent已创建")
    return app