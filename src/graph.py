from langgraph.graph import StateGraph, START, END
from .state import AgentState
from .node.normal_node import normal_node
from .node.user_communication_node import user_communication_node
from .node.plan_node import plan_node
from .node.detect_node import detect_node
from .node.prepare_node import prepare_node
from .node.user_confirm_node import user_confirm_node
from .node.execute_node import execute_node
from .route import route_from_user_communication, route_after_confirm

def create_agent():
    """创建LangGraph Agent - User Communication为主的Plan-and-Execute架构"""
    
    # 创建状态图
    workflow = StateGraph(AgentState)
    
    # 添加所有节点
    workflow.add_node("normal", normal_node)
    workflow.add_node("user_communication", user_communication_node)
    workflow.add_node("plan", plan_node)
    workflow.add_node("detect", detect_node)
    workflow.add_node("prepare", prepare_node)
    workflow.add_node("user_confirm", user_confirm_node)
    workflow.add_node("execute", execute_node)
    
    # 入口点：直接进入User Communication节点
    workflow.add_edge(START, "user_communication")
    
    # User Communication节点的条件路由
    workflow.add_conditional_edges(
        "user_communication", 
        route_from_user_communication,
        {
            "end": END,             # 结束流程
            "normal": "normal",     # 进入意图分析
            "plan": "plan"          # 进入分析流程
        }
    )
    
    # Normal节点路由（仅回到user_communication）
    workflow.add_edge("normal", "user_communication")
    
    # 分析流程: Plan → Detect → Prepare (简化路由)
    workflow.add_edge("plan", "detect")
    workflow.add_edge("detect", "prepare")  # Detect完成所有任务后直接进入Prepare
    workflow.add_edge("prepare", "user_confirm")
    
    # 用户确认后的路由
    workflow.add_conditional_edges(
        "user_confirm",
        route_after_confirm,
        {
            "execute": "execute",
            "plan": "plan",  # replan直接回到plan节点重新规划
            "cancel": "user_communication",
            "quit": END
        }
    )
    
    # 执行完成后结束
    workflow.add_edge("execute", END)
    
    # 编译图
    app = workflow.compile()
    
    print("🤖 RNA-seq智能分析助手已启动")
    print("   架构: User Communication → Normal → Plan → Execute (支持循环重规划)")
    return app