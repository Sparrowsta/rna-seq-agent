"""
路由决策器 - 基于执行结果推断下一步行动

统一规则：
- 失败: 返回用户确认
- 成功: 根据执行模式决定是否继续
"""

from .state import AgentState
from .logging_bootstrap import get_logger

logger = get_logger("rna.route_decider")


def compute_action(mode: str, step: str, success: bool) -> str:
    """
    统一路由决策函数
    
    Args:
        mode: 执行模式 (single/optimized/batch_optimize/yolo)
        step: 当前步骤 (fastp/star/hisat2/featurecounts)  
        success: 执行是否成功
        
    Returns:
        决策结果: continue_next/return_confirm
    """
    mode = (mode or '').lower()
    
    logger.debug(f"[COMPUTE-ACTION] 策略决策: mode={mode}, step={step}, success={success}")
    
    # 失败优先：失败统一返回确认
    if not success:
        logger.info(f"[COMPUTE-ACTION] 失败处理: {step} -> return_confirm")
        return 'return_confirm'
    
    # optimized：每步都确认（即使成功也确认）
    if mode == 'optimized':
        logger.info(f"[COMPUTE-ACTION] OPTIMIZED模式: {step} 成功 -> return_confirm")
        return 'return_confirm'
    
    # batch：末端（featurecounts）确认，其余自动推进
    if mode == 'batch_optimize' and step == 'featurecounts':
        logger.info(f"[COMPUTE-ACTION] BATCH_OPTIMIZE模式: {step} 末端确认 -> return_confirm")
        return 'return_confirm'
    
    # single/yolo/batch(非末端)：成功即前进
    logger.info(f"[COMPUTE-ACTION] {mode.upper()}模式: {step} 成功 -> continue_next")
    return 'continue_next'


def decide_next_action_fastp(state: AgentState) -> str:
    """
    决策FastP后的下一步行动
    
    Returns:
        - "continue_next": 继续比对步骤
        - "return_confirm": 返回用户确认
    """
    # 提取状态信息
    fastp_results = getattr(state, 'fastp_results', {}) or {}
    
    success = bool(fastp_results.get('success', False))
    
    logger.debug(f"[ROUTE-DECIDER] FastP - mode: {getattr(state, 'execution_mode', 'single')}, success: {success}")
    
    # 调用极简的统一策略函数
    return compute_action(
        mode=getattr(state, 'execution_mode', 'single'),
        step='fastp',
        success=success
    )


def decide_next_action_star(state: AgentState) -> str:
    """
    决策STAR后的下一步行动
    
    Returns:
        - "continue_next": 继续定量步骤
        - "return_confirm": 返回用户确认
    """
    # 提取状态信息
    star_results = getattr(state, 'star_results', {}) or {}
    
    success = bool(star_results.get('success', False))
    
    logger.debug(f"[ROUTE-DECIDER] STAR - mode: {getattr(state, 'execution_mode', 'single')}, success: {success}")
    
    # 调用极简的统一策略函数
    return compute_action(
        mode=getattr(state, 'execution_mode', 'single'),
        step='star',
        success=success
    )


def decide_next_action_hisat2(state: AgentState) -> str:
    """
    决策HISAT2后的下一步行动
    
    Returns:
        - "continue_next": 继续定量步骤
        - "return_confirm": 返回用户确认
    """
    # 提取状态信息
    hisat2_results = getattr(state, 'hisat2_results', {}) or {}
    
    success = bool(hisat2_results.get('success', False))
    
    logger.debug(f"[ROUTE-DECIDER] HISAT2 - mode: {getattr(state, 'execution_mode', 'single')}, success: {success}")
    
    # 调用极简的统一策略函数
    return compute_action(
        mode=getattr(state, 'execution_mode', 'single'),
        step='hisat2',
        success=success
    )


def decide_next_action_featurecounts(state: AgentState) -> str:
    """
    决策FeatureCounts后的下一步行动
    
    Returns:
        - "continue_next": 继续分析步骤
        - "return_confirm": 返回用户确认  
    """
    # 提取状态信息
    featurecounts_results = getattr(state, 'featurecounts_results', {}) or {}
    
    success = bool(featurecounts_results.get('success', False))
    
    logger.debug(f"[ROUTE-DECIDER] FeatureCounts - mode: {getattr(state, 'execution_mode', 'single')}, success: {success}")
    
    # 调用极简的统一策略函数
    return compute_action(
        mode=getattr(state, 'execution_mode', 'single'),
        step='featurecounts',
        success=success
    )
