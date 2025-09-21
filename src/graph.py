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
from .node.hisat2_node import hisat2_node
from .node.featurecounts_node import featurecounts_node
from .node.analysis_node import analysis_node
from .route import (
    route_from_user_communication, 
    route_after_confirm, 
    route_after_fastp,
    route_after_star,
    route_after_hisat2,
    route_after_featurecount,
    route_after_analysis
)
from .logging_bootstrap import get_logger, is_debug_enabled

logger = get_logger("rna.graph")

def create_agent():
    """创建LangGraph Agent - 支持STAR/HISAT2双比对器的完整RNA-seq流程"""
    
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
    workflow.add_node("hisat2", hisat2_node)
    workflow.add_node("featurecounts", featurecounts_node)
    workflow.add_node("analysis", analysis_node)
    workflow.add_node("modify", modify_node)
    
    # 入口点：直接进入User Communication节点
    workflow.add_edge(START, "user_communication")
    
    # User Communication节点的条件路由
    workflow.add_conditional_edges(
        "user_communication", 
        route_from_user_communication,
        {
            "end": END,                     # 结束流程
            "normal": "normal",             # 进入意图分析
            "detect": "detect",             # 进入检测流程（去除Plan节点）
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
            "fastp": "fastp",                     # 开始FastP处理
            "star": "star",                       # 继续STAR比对
            "hisat2": "hisat2",                   # 继续HISAT2比对
            "featurecounts": "featurecounts",     # 继续FeatureCounts定量
            "analysis": "analysis",               # 继续综合分析
            "continue_star": "star",              # /continue命令：继续到STAR
            "continue_hisat2": "hisat2",          # /continue命令：继续到HISAT2
            "continue_featurecounts": "featurecounts",  # /continue命令：继续到FeatureCounts
            "continue_analysis": "analysis",      # /continue命令：继续到Analysis
            "modify": "modify",                   # 修改配置路由
            "user_confirm": "user_confirm",       # 回退到确认（兜底，防未识别输入）
            "cancel": "user_communication",
            "quit": END
        }
    )
    
    # Modify节点完成后直接返回User Confirm节点
    workflow.add_edge("modify", "user_confirm")
    
    # FastP节点完成后的路由：根据配置和mode决定进入STAR或HISAT2
    workflow.add_conditional_edges(
        "fastp",
        route_after_fastp,
        {
            "star": "star",                       # 进入STAR比对
            "hisat2": "hisat2",                   # 进入HISAT2比对
            "user_confirm": "user_confirm",       # 回到确认（优化模式）
        }
    )
    
    # STAR节点完成后的路由
    workflow.add_conditional_edges(
        "star",
        route_after_star,
        {
            "featurecounts": "featurecounts",   # 继续FeatureCount定量
            "user_confirm": "user_confirm",   # 回到确认（优化模式或错误） 
        }
    )
    
    # HISAT2节点完成后的路由
    workflow.add_conditional_edges(
        "hisat2",
        route_after_hisat2,
        {
            "featurecounts": "featurecounts",   # 继续FeatureCount定量
            "user_confirm": "user_confirm",   # 回到确认（优化模式或错误）
        }
    )
    
    # FeatureCount节点完成后的路由
    workflow.add_conditional_edges(
        "featurecounts", 
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
            "user_confirm": "user_confirm",  # 返回用户确认界面
        }
    )
    
    # 编译图
    app = workflow.compile()
    
    # 如果启用调试模式，记录调试信息
    if is_debug_enabled():
        logger.debug("调试模式已启用，支持LangGraph stream调试")
        
    return app