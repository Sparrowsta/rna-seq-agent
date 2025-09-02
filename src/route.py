from langgraph.graph import END
from .state import AgentState

def route_from_user_communication(state: AgentState) -> str:
    """User Communication节点后的路由决策"""
    routing_decision = state.routing_decision
    
    if routing_decision == "plan":
        print("🚀 进入Plan分析流程")
        return "plan"
    elif routing_decision == "normal":
        print("🧠 进入意图分析")
        return "normal"
    elif routing_decision == "end":
        print("🔚 会话结束")
        return "end"
    else:
        print(f"⚠️ 未知路由决策: {routing_decision}，默认结束会话")
        return "end"

def should_continue(state: AgentState) -> str:
    """决定是否继续执行"""
    plan = state.plan
    # plan现在是List[List[str]]格式，检查是否有任务组
    if plan and any(group for group in plan):
        return "detect"
    else:
        return "prepare"

def route_after_confirm(state: AgentState) -> str:
    """用户确认后的路由决策"""
    user_decision = state.user_decision.lower() if state.user_decision else ""
    
    print(f"\n🔍 [DEBUG] 路由决策分析:")
    print(f"   用户决策: '{state.user_decision}'")
    print(f"   标准化后: '{user_decision}'")
    
    if user_decision == "execute":
        print("🚀 [ROUTE] 用户选择执行分析")
        return "execute"
    elif user_decision in ["modify", "replan"]:
        print("🔄 [ROUTE] 用户选择重新规划，回到plan节点") 
        return "plan"
    elif user_decision == "cancel":
        print("❌ [ROUTE] 用户选择取消分析")
        return "cancel"
    elif user_decision == "quit":
        print("🚪 [ROUTE] 用户选择退出程序")
        return "quit"
    else:
        print(f"⚠️ [ROUTE] 未识别的决策 '{user_decision}'，默认取消")
        return "cancel"

def route_after_analysis(state: AgentState) -> str:
    """Analysis节点分析完毕后的路由决策"""
    print("✅ [ROUTE] 分析总结完成，返回用户交互")
    return "user_communication"