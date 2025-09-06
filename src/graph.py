from langgraph.graph import StateGraph, START, END
from .state import AgentState
from .node.normal_node import normal_node
from .node.user_communication_node import user_communication_node
from .node.detect_node import detect_node
from .node.prepare_node import prepare_node
from .node.user_confirm_node import user_confirm_node
from .node.fastp_node import fastp_node
from .node.modify_node import modify_node
from .route import route_from_user_communication, route_after_confirm, route_after_fastp

def create_agent():
    """创建LangGraph Agent - User Communication为主的Plan-and-Execute架构"""
    
    # 创建状态图
    workflow = StateGraph(AgentState)
    
    # 添加所有节点
    workflow.add_node("normal", normal_node)
    workflow.add_node("user_communication", user_communication_node)
    workflow.add_node("detect", detect_node)
    workflow.add_node("prepare", prepare_node)
    workflow.add_node("user_confirm", user_confirm_node)
    workflow.add_node("fastp", fastp_node)
    workflow.add_node("modify", modify_node)  # 添加modify节点
    
    # 入口点：直接进入User Communication节点
    workflow.add_edge(START, "user_communication")
    
    # User Communication节点的条件路由
    workflow.add_conditional_edges(
        "user_communication", 
        route_from_user_communication,
        {
            "end": END,             # 结束流程
            "normal": "normal",     # 进入意图分析
            "detect": "detect"         # 进入检测流程（去除Plan节点）
        }
    )
    
    # Normal节点路由（仅回到user_communication）
    workflow.add_edge("normal", "user_communication")
    
    # 分析流程: 直接 Detect → Prepare
    workflow.add_edge("detect", "prepare")  # Detect完成所有任务后直接进入Prepare
    workflow.add_edge("prepare", "user_confirm")
    
    # 用户确认后的路由
    workflow.add_conditional_edges(
        "user_confirm",
        route_after_confirm,
        {
            "fastp": "fastp",                 # 统一执行路由：所有分析任务都通过fastp处理
            "modify": "modify",               # 修改配置路由 - 先进入modify节点
            "cancel": "user_communication",
            "quit": END
        }
    )
    
    # Modify节点完成后直接返回User Confirm节点
    workflow.add_edge("modify", "user_confirm")
    
    # FastP节点完成后：单次执行直接结束；优化执行回到确认
    workflow.add_conditional_edges(
        "fastp",
        route_after_fastp,
        {
            "user_confirm": "user_confirm",
            "end": END
        }
    )
    
    # 编译图
    app = workflow.compile()
    
    print("🤖 RNA-seq智能分析助手已启动")
    print("   架构: User Communication → Normal → Detect → Prepare → Confirm → (Modify →) FastP → (END/Confirm)")
    return app
