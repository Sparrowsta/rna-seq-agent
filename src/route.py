from langgraph.graph import END
from .state import AgentState
from .logging_bootstrap import get_logger
from .route_decider import (
    decide_next_action_fastp,
    decide_next_action_star,
    decide_next_action_hisat2,
    decide_next_action_featurecounts
)

logger = get_logger("rna.route")

def route_from_user_communication(state: AgentState) -> str:
    """User Communication节点后的路由决策"""
    routing_decision = state.routing_decision
    
    if routing_decision == "plan":
        logger.info("进入检测流程")
        return "detect"
    elif routing_decision == "execute":  # 新增execute模式支持
        logger.info("进入执行模式检测流程")
        return "detect"
    elif routing_decision == "normal":
        logger.info("进入意图分析")
        return "normal"
    elif routing_decision == "end":
        logger.info("会话结束")
        return "end"
    else:
        logger.warning(f"未知路由决策: {routing_decision}，请重新输入")
        return "normal"


def route_after_confirm(state: AgentState) -> str:
    """用户确认后的路由决策"""
    user_decision = state.user_decision.lower() if state.user_decision else ""
    
    logger.debug(f"[DEBUG] 路由决策分析:")
    logger.debug(f"用户决策: '{state.user_decision}'")
    logger.debug(f"标准化后: '{user_decision}'")
    
    if user_decision == "execute":
        logger.info("[ROUTE] 用户选择执行分析")
        logger.info("[ROUTE] 统一路由到FastP节点处理")
        return "fastp"
    elif user_decision in {"fastp", "star", "featurecounts"}:
        # 直接回到指定步骤对应的Agent（用于 /re_opt 二次优化等）
        logger.info(f"[ROUTE] 返回到当前步骤进行二次优化: {user_decision}")
        return user_decision
    elif user_decision == "continue_star":
        logger.info("[ROUTE] 继续到STAR比对")
        return "star"
    elif user_decision == "continue_hisat2":
        logger.info("[ROUTE] 继续到HISAT2比对")
        return "hisat2"
    elif user_decision == "continue_featurecounts":
        logger.info("[ROUTE] 继续到FeatureCounts定量")
        return "featurecounts"
    elif user_decision == "continue_analysis":
        logger.info("[ROUTE] 继续到综合分析")
        return "analysis"
    elif user_decision == "modify":
        logger.info("[ROUTE] 用户选择修改配置")
        return "modify"
    elif user_decision == "cancel":
        logger.info("[ROUTE] 用户选择取消分析")
        return "cancel"
    elif user_decision == "quit":
        logger.info("[ROUTE] 用户选择退出程序")
        return "quit"
    else:
        logger.warning(f"[ROUTE] 未识别的决策 '{user_decision}'，请重新选择")
        return "user_confirm"

def route_after_fastp(state: AgentState) -> str:
    """FastP节点后的路由：
    根据决策器返回的 next_action 选择下一步：
    - continue_next: 根据配置选择比对器（star/hisat2）继续执行
    - return_confirm: 返回用户确认界面（失败或有优化建议）
    - re_run: 重新运行FastP（可自动修复的错误）
    """
    # 调用路由决策器获取下一步行动
    next_action = decide_next_action_fastp(state)
    
    logger.debug(f"[ROUTE-FASTP] 决策器返回: {next_action}")
    
    if next_action == "return_confirm":
        logger.info("[ROUTE] FastP需要用户确认（失败或有优化建议），返回确认界面")
        return "user_confirm"
    
    elif next_action == "re_run":
        logger.info("[ROUTE] FastP检测到可自动修复错误，重新运行")
        return "fastp"
    
    elif next_action == "continue_next":
        # 选择比对器 - 优先级：配置 > 用户需求 > 默认(star)
        # 支持两种键名：align_tool（首选）与 aligner（兼容旧字段）
        aligner = "star"  # 默认使用STAR
        if hasattr(state, 'nextflow_config') and state.nextflow_config:
            cfg = state.nextflow_config
            aligner = (cfg.get('align_tool')
                       or cfg.get('aligner')
                       or 'star')
        elif hasattr(state, 'user_requirements') and state.user_requirements:
            req = state.user_requirements
            aligner = (req.get('align_tool')
                       or req.get('aligner')
                       or 'star')
        
        # 确保比对器名称标准化
        aligner = aligner.lower()
        if aligner not in ['star', 'hisat2']:
            logger.warning(f"[ROUTE] 未知比对器 '{aligner}'，使用默认STAR")
            aligner = 'star'

        logger.info(f"[ROUTE] FastP完成，继续{aligner.upper()}比对")
        return aligner
    
    else:
        # 兜底逻辑：未知决策返回确认界面
        logger.warning(f"[ROUTE] 未知决策 '{next_action}'，返回确认界面")
        return "user_confirm"


def route_after_star(state: AgentState) -> str:
    """STAR节点后的路由：
    根据决策器返回的 next_action 选择下一步：
    - continue_next: 继续FeatureCounts定量
    - return_confirm: 返回用户确认界面（失败或有优化建议）
    - re_run: 重新运行STAR（可自动修复的错误）
    """
    # 调用路由决策器获取下一步行动
    next_action = decide_next_action_star(state)
    
    logger.debug(f"[ROUTE-STAR] 决策器返回: {next_action}")
    
    if next_action == "return_confirm":
        logger.info("[ROUTE] STAR需要用户确认（失败或有优化建议），返回确认界面")
        return "user_confirm"
    
    elif next_action == "re_run":
        logger.info("[ROUTE] STAR检测到可自动修复错误，重新运行")
        return "star"
    
    elif next_action == "continue_next":
        logger.info("[ROUTE] STAR比对成功，继续FeatureCounts定量")
        return "featurecounts"
    
    else:
        # 兜底逻辑：未知决策返回确认界面
        logger.warning(f"[ROUTE] 未知决策 '{next_action}'，返回确认界面")
        return "user_confirm"
       

def route_after_hisat2(state: AgentState) -> str:
    """HISAT2节点后的路由：
    根据决策器返回的 next_action 选择下一步：
    - continue_next: 继续FeatureCounts定量
    - return_confirm: 返回用户确认界面（失败或有优化建议）
    - re_run: 重新运行HISAT2（可自动修复的错误）
    """
    # 调用路由决策器获取下一步行动
    next_action = decide_next_action_hisat2(state)
    
    logger.debug(f"[ROUTE-HISAT2] 决策器返回: {next_action}")
    
    if next_action == "return_confirm":
        logger.info("[ROUTE] HISAT2需要用户确认（失败或有优化建议），返回确认界面")
        return "user_confirm"
    
    elif next_action == "re_run":
        logger.info("[ROUTE] HISAT2检测到可自动修复错误，重新运行")
        return "hisat2"
    
    elif next_action == "continue_next":
        logger.info("[ROUTE] HISAT2比对成功，继续FeatureCounts定量")
        return "featurecounts"
    
    else:
        # 兜底逻辑：未知决策返回确认界面
        logger.warning(f"[ROUTE] 未知决策 '{next_action}'，返回确认界面")
        return "user_confirm"



def route_after_featurecount(state: AgentState) -> str:
    """FeatureCount节点后的路由：
    根据决策器返回的 next_action 选择下一步：
    - continue_next: 进入综合分析
    - return_confirm: 返回用户确认界面（失败或有优化建议）
    - re_run: 重新运行FeatureCounts（可自动修复的错误）
    """
    # 调用路由决策器获取下一步行动
    next_action = decide_next_action_featurecounts(state)
    
    logger.debug(f"[ROUTE-FEATURECOUNTS] 决策器返回: {next_action}")
    
    if next_action == "return_confirm":
        logger.info("[ROUTE] FeatureCounts需要用户确认（失败或有优化建议），返回确认界面")
        return "user_confirm"
    
    elif next_action == "re_run":
        logger.info("[ROUTE] FeatureCounts检测到可自动修复错误，重新运行")
        return "featurecounts"
    
    elif next_action == "continue_next":
        logger.info("[ROUTE] FeatureCounts定量成功，进入综合分析")
        return "analysis"
    
    else:
        # 兜底逻辑：未知决策返回确认界面
        logger.warning(f"[ROUTE] 未知决策 '{next_action}'，返回确认界面")
        return "user_confirm"

    
def route_after_analysis(state: AgentState) -> str:
    """Analysis节点后的路由 - 基于四模式路由策略v1：
    - single/yolo模式：直接结束（END）
    - optimized/batch_optimize模式：返回用户确认界面
    """
    analysis_success = getattr(state, 'success', False)
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()

    if not analysis_success:
        logger.info(f"[ROUTE] 分析失败（模式：{mode}），返回用户确认界面")
        return "user_confirm"

    if mode in ('single', 'yolo'):
        logger.info(f"[ROUTE] {mode.upper()}模式：分析完成，流程结束")
        return END
    if mode in ('optimized', 'batch_optimize'):
        logger.info(f"[ROUTE] {mode.upper()}模式：分析完成，返回用户确认界面")
        return "user_confirm"

    logger.warning(f"[ROUTE] 未知执行模式 '{mode}'，兜底返回确认界面")
    return "user_confirm"
