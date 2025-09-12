from langgraph.graph import END
from .state import AgentState

def route_from_user_communication(state: AgentState) -> str:
    """User Communication节点后的路由决策"""
    routing_decision = state.routing_decision
    
    if routing_decision == "plan":
        print("🚀 进入检测流程")
        return "detect"
    elif routing_decision == "normal":
        print("🧠 进入意图分析")
        return "normal"
    elif routing_decision == "end":
        print("🔚 会话结束")
        return "end"
    else:
        print(f"⚠️ 未知路由决策: {routing_decision}，请重新输入")
        return "normal"


def route_after_confirm(state: AgentState) -> str:
    """用户确认后的路由决策"""
    user_decision = state.user_decision.lower() if state.user_decision else ""
    
    print(f"\n🔍 [DEBUG] 路由决策分析:")
    print(f"   用户决策: '{state.user_decision}'")
    print(f"   标准化后: '{user_decision}'")
    
    if user_decision == "execute":
        print("🚀 [ROUTE] 用户选择执行分析")
        print("🧬 [ROUTE] 统一路由到FastP节点处理")
        return "fastp"
    elif user_decision in {"fastp", "star", "featurecounts"}:
        # 直接回到指定步骤对应的Agent（用于 /re_opt 二次优化等）
        print(f"♻️ [ROUTE] 返回到当前步骤进行二次优化: {user_decision}")
        return user_decision
    elif user_decision == "continue_star":
        print("🎯 [ROUTE] 继续到STAR比对")
        return "star"
    elif user_decision == "continue_featurecounts":
        print("📊 [ROUTE] 继续到FeatureCounts定量")
        return "featurecounts"
    elif user_decision == "continue_analysis":
        print("📈 [ROUTE] 继续到综合分析")
        return "analysis"
    elif user_decision == "modify":
        print("🔧 [ROUTE] 用户选择修改配置")
        return "modify"
    elif user_decision == "cancel":
        print("❌ [ROUTE] 用户选择取消分析")
        return "cancel"
    elif user_decision == "quit":
        print("🚪 [ROUTE] 用户选择退出程序")
        return "quit"
    else:
        print(f"⚠️ [ROUTE] 未识别的决策 '{user_decision}'，请重新选择")
        return "user_confirm"

def route_after_fastp(state: AgentState) -> str:
    """FastP节点后的路由：
    - 单次执行（single）：继续STAR比对
    - 优化执行（optimized）：回到用户确认进行参数微调
    - 批次优化（batch_optimize）：继续STAR比对（收集优化建议但不中断）
    - 错误情况：回到用户确认
    """
    mode = (getattr(state, 'execution_mode', 'single') or 'single').lower()
    fastp_results = getattr(state, 'fastp_results', {}) or {}

    fastp_success = fastp_results.get('success', False)
    print(f"\n🔎 [ROUTE-FASTP] mode={mode}, success={fastp_success}")
    if not fastp_success:
        print("❌ [ROUTE] FastP未成功或结果缺失，返回确认界面")
        print(f"   [DEBUG] fastp_results keys: {list(fastp_results.keys())}")
        return "user_confirm"

    if mode == 'yolo':
        print("🎯 [ROUTE] YOLO模式：FastP完成后自动进入STAR比对")
        return "star"
    elif mode == 'optimized':
        print("🔁 [ROUTE] 优化执行模式：FastP 完成后返回确认进行参数微调")
        return "user_confirm"
    elif mode == 'batch_optimize':
        print("📊 [ROUTE] 批次优化模式：FastP 完成后继续STAR比对（收集优化建议）")
        return "star"
    else:  # single 及其他
        print("🧬 [ROUTE] 单次执行模式：FastP 完成后继续STAR比对")
        return "star"


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
    print(f"\n🔎 [ROUTE-STAR] mode={mode}, success={star_success}")
    
    # 检查STAR是否成功完成
    if star_success:
        if mode == 'yolo':
            print("🎯 [ROUTE] YOLO模式：STAR完成后自动进入FeatureCount定量")
            return "featurecounts"
        elif mode == 'optimized':
            print("🔁 [ROUTE] 优化执行模式：STAR完成后返回确认进行参数微调")
            return "user_confirm"
        else:  # single 或 batch_optimize 等
            print("🧬 [ROUTE] STAR比对成功，继续FeatureCount定量")
            return "featurecounts"
    else:
        print("❌ [ROUTE] star执行失败，返回确认界面")
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
    print(f"\n🔎 [ROUTE-FEATURECOUNTS] mode={mode}, success={fc_success}")

    # 检查FeatureCounts是否成功完成
    if fc_success:
        if mode == 'yolo':
            print("🎯 [ROUTE] YOLO模式：FeatureCount完成后自动进入综合分析")
            return "analysis"
        elif mode == 'optimized':
            print("🔁 [ROUTE] 优化执行模式：FeatureCount完成后返回确认进行参数微调")
            return "user_confirm"
        elif mode == 'batch_optimize':
            print("📊 [ROUTE] 批量优化模式：FeatureCount完成，返回确认界面显示优化建议")
            return "user_confirm"
        else:  # single 或其他
            print("🧬 [ROUTE] FeatureCount定量成功，进入综合分析")
            return "analysis"
    else:
        print("❌ [ROUTE] FeatureCount定量失败，返回确认界面")
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
            print("🧬 [ROUTE] 满足分析条件，进入综合分析")
            return "analysis"
    
    print("❌ [ROUTE] 不满足分析条件，返回确认界面")
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
            print("🎯 [ROUTE] YOLO模式：分析完成，流程结束")
            return END
        else:
            print("✅ [ROUTE] 分析完成，返回用户确认界面")
            return "user_confirm"
    else:
        print("❌ [ROUTE] 分析失败，返回用户确认界面")
        return "user_confirm"
