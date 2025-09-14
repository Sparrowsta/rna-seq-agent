from langgraph.graph import END
from .state import AgentState
from .logging_bootstrap import get_logger

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
    - 检查配置中的比对器选择（align_tool/aligner: star 或 hisat2）
    - 单次执行（single）：根据配置进入相应比对
    - 优化执行（optimized）：回到用户确认进行参数微调
    - 批次优化（batch_optimize）：根据配置进入相应比对（收集优化建议但不中断）
    - 错误情况：回到用户确认
    """
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    fastp_results = getattr(state, 'fastp_results', {}) or {}

    fastp_success = fastp_results.get('success', False)
    logger.debug(f"[ROUTE-FASTP] mode={mode}, success={fastp_success}")
    if not fastp_success:
        logger.info("[ROUTE] FastP未成功或结果缺失，返回确认界面")
        logger.debug(f"[DEBUG] fastp_results keys: {list(fastp_results.keys())}")
        return "user_confirm"

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

    logger.info(f"[ROUTE] 选择的比对器: {aligner.upper()}")

    if mode == 'yolo':
        logger.info(f"[ROUTE] YOLO模式：FastP完成后自动进入{aligner.upper()}比对")
        return aligner
    elif mode == 'optimized':
        logger.info("[ROUTE] 优化执行模式：FastP 完成后返回确认进行参数微调")
        return "user_confirm"
    elif mode == 'batch_optimize':
        logger.info(f"[ROUTE] 批次优化模式：FastP 完成后继续{aligner.upper()}比对（收集优化建议）")
        return aligner
    else:  # single 及其他
        logger.info(f"[ROUTE] 单次执行模式：FastP 完成后继续{aligner.upper()}比对")
        return aligner


def route_after_star(state: AgentState) -> str:
    """STAR节点后的路由：
    - 单次执行：继续FeatureCount定量
    - 优化执行：回到用户确认进行参数微调
    - 批次优化模式：继续FeatureCount（收集优化建议但不中断）
    - 其他错误：回到用户确认
    """
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    star_results = getattr(state, 'star_results', {}) or {}

    star_success = star_results.get('success', False)
    logger.debug(f"[ROUTE-STAR] mode={mode}, success={star_success}")
    
    # 检查STAR是否成功完成
    if star_success:
        if mode == 'yolo':
            logger.info("[ROUTE] YOLO模式：STAR完成后自动进入FeatureCount定量")
            return "featurecounts"
        elif mode == 'optimized':
            logger.info("[ROUTE] 优化执行模式：STAR完成后返回确认进行参数微调")
            return "user_confirm"
        else:  # single 或 batch_optimize 等
            logger.info("[ROUTE] STAR比对成功，继续FeatureCount定量")
            return "featurecounts"
    else:
        logger.info("[ROUTE] star执行失败，返回确认界面")
        return "user_confirm"
       

def route_after_hisat2(state: AgentState) -> str:
    """HISAT2节点后的路由：
    - 单次执行：继续FeatureCount定量
    - 优化执行：回到用户确认进行参数微调
    - 批次优化模式：继续FeatureCount（收集优化建议但不中断）
    - 其他错误：回到用户确认
    """
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    hisat2_results = getattr(state, 'hisat2_results', {}) or {}

    hisat2_success = hisat2_results.get('success', False)
    logger.debug(f"[ROUTE-HISAT2] mode={mode}, success={hisat2_success}")
    
    # 检查HISAT2是否成功完成
    if hisat2_success:
        if mode == 'yolo':
            logger.info("[ROUTE] YOLO模式：HISAT2完成后自动进入FeatureCount定量")
            return "featurecounts"
        elif mode == 'optimized':
            logger.info("[ROUTE] 优化执行模式：HISAT2完成后返回确认进行参数微调")
            return "user_confirm"
        else:  # single 或 batch_optimize 等
            logger.info("[ROUTE] HISAT2比对成功，继续FeatureCount定量")
            return "featurecounts"
    else:
        logger.info("[ROUTE] hisat2执行失败，返回确认界面")
        return "user_confirm"



def route_after_featurecount(state: AgentState) -> str:
    """FeatureCount节点后的路由：
    - 成功：进入综合分析
    - 优化模式错误：回到用户确认
    - 批量优化模式：继续分析（收集所有优化建议）
    - 其他错误：回到用户确认
    """
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    featurecounts_results = getattr(state, 'featurecounts_results', {}) or {}

    fc_success = featurecounts_results.get('success', False)
    logger.debug(f"[ROUTE-FEATURECOUNTS] mode={mode}, success={fc_success}")

    # 检查FeatureCounts是否成功完成
    if fc_success:
        if mode == 'yolo':
            logger.info("[ROUTE] YOLO模式：FeatureCount完成后自动进入综合分析")
            return "analysis"
        elif mode == 'optimized':
            logger.info("[ROUTE] 优化执行模式：FeatureCount完成后返回确认进行参数微调")
            return "user_confirm"
        elif mode == 'batch_optimize':
            logger.info("[ROUTE] 批量优化模式：FeatureCount完成，返回确认界面显示优化建议")
            return "user_confirm"
        else:  # single 或其他
            logger.info("[ROUTE] FeatureCount定量成功，进入综合分析")
            return "analysis"
    else:
        logger.info("[ROUTE] FeatureCount定量失败，返回确认界面")
        return "user_confirm"


def route_to_analysis(state: AgentState) -> str:
    """检查是否可以进入Analysis节点：
    - 满足条件：进入分析节点
    - 不满足：回到用户确认
    """
    mode = getattr(state, 'execution_mode', '')
    
    # 检查是否有分析执行请求
    if mode == "single_run":
        # 检查是否有必要的数据进行分析
        fastp_results = getattr(state, 'fastp_results', {})
        star_results = getattr(state, 'star_results', {}) 
        featurecounts_results = getattr(state, 'featurecounts_results', {})
        
        if fastp_results or star_results or featurecounts_results:
            logger.info("[ROUTE] 满足分析条件，进入综合分析")
            return "analysis"
    
    logger.info("[ROUTE] 不满足分析条件，返回确认界面")
    return "user_confirm"


def route_after_analysis(state: AgentState) -> str:
    """Analysis节点后的路由：
    - YOLO模式：直接结束
    - 其他模式：返回用户确认界面
    """
    # 读取节点顶层success字段，符合success-first约定
    analysis_success = getattr(state, 'success', False)
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    
    if analysis_success:
        if mode == 'yolo':
            logger.info("[ROUTE] YOLO模式：分析完成，流程结束")
            return END
        else:
            logger.info("[ROUTE] 分析完成，返回用户确认界面")
            return "user_confirm"
    else:
        logger.info("[ROUTE] 分析失败，返回用户确认界面")
        return "user_confirm"
