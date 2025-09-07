from langgraph.graph import StateGraph, START, END
from .state import AgentState
from .node.normal_node import normal_node
from .node.user_communication_node import user_communication_node
from .node.detect_node import detect_node
from .node.prepare_node import prepare_node
from .node.user_confirm_node import user_confirm_node
from .node.fastp_node import fastp_node
from .node.modify_node import modify_node
from .node.star_node import star_node
from .node.featurecount_node import featurecount_node
from .node.analysis_node import analysis_node
from .route import (
    route_from_user_communication, 
    route_after_confirm, 
    route_after_fastp,
    route_after_star,
    route_after_featurecount,
    route_to_analysis,
    route_after_analysis
)

def create_agent():
    """创建LangGraph Agent - 支持STAR-FeatureCount-Analysis完整流程"""
    
    # 创建状态图
    workflow = StateGraph(AgentState)
    
    # 添加所有节点
    workflow.add_node("normal", normal_node)
    workflow.add_node("user_communication", user_communication_node)
    workflow.add_node("detect", detect_node)
    workflow.add_node("prepare", prepare_node)
    workflow.add_node("user_confirm", user_confirm_node)
    workflow.add_node("fastp", fastp_node)
    workflow.add_node("star", star_node)
    workflow.add_node("featurecount", featurecount_node)
    workflow.add_node("analysis", analysis_node)
    workflow.add_node("modify", modify_node)
    
    # 入口点：直接进入User Communication节点
    workflow.add_edge(START, "user_communication")
    
    # User Communication节点的条件路由
    workflow.add_conditional_edges(
        "user_communication", 
        route_from_user_communication,
        {
            "end": END,             # 结束流程
            "normal": "normal",     # 进入意图分析
            "detect": "detect"      # 进入检测流程（去除Plan节点）
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
            "fastp": "fastp",                 # 开始FastP处理
            "modify": "modify",               # 修改配置路由
            "cancel": "user_communication",
            "quit": END
        }
    )
    
    # Modify节点完成后直接返回User Confirm节点
    workflow.add_edge("modify", "user_confirm")
    
    # FastP节点完成后的路由：根据mode决定下一步
    workflow.add_conditional_edges(
        "fastp",
        route_after_fastp,
        {
            "star": "star",                   # 继续STAR比对
            "user_confirm": "user_confirm",   # 回到确认（优化模式）
        }
    )
    
    # STAR节点完成后的路由
    workflow.add_conditional_edges(
        "star",
        route_after_star,
        {
            "featurecount": "featurecount",   # 继续FeatureCount定量
            "user_confirm": "user_confirm",   # 回到确认（优化模式或错误） 
        }
    )
    
    # FeatureCount节点完成后的路由
    workflow.add_conditional_edges(
        "featurecount", 
        route_after_featurecount,
        {
            "analysis": "analysis",           # 进入综合分析
            "user_confirm": "user_confirm",   # 回到确认（优化模式或错误）
        }
    )
    
    # Analysis节点完成后的路由
    workflow.add_conditional_edges(
        "analysis",
        route_after_analysis,
        {
            "user_communication": "user_communication",  # 返回用户交互
        }
    )
    
    # 编译图
    app = workflow.compile()
    
    print("🤖 RNA-seq智能分析助手已启动")
    print("   架构: User Communication → Normal → Detect → Prepare → Confirm")
    print("   流程: (Modify →) FastP → STAR → FeatureCount → Analysis → (END/Confirm)")
    return app
