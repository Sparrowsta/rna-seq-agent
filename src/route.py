from langgraph.graph import END
from .state import UserCommunicationNodeState, DetectNodeState, UserConfirmState

async def route_from_user_communication(state: UserCommunicationNodeState) -> str:
    """User Communication节点后的路由决策"""
    routing_decision = state.get("routing_decision", "end")
    
    if routing_decision == "plan":
        print("🚀 进入Plan分析流程")
        return "plan"
    elif routing_decision == "normal":
        print("🧠 进入意图分析")
        return "normal"
    else:
        print("🔚 会话结束")
        return "end"

async def should_continue(state: DetectNodeState) -> str:
    """决定是否继续执行"""
    plan = state.get("plan", [])
    if plan:
        return "detect"
    else:
        return "prepare"

async def route_after_confirm(state: UserConfirmState) -> str:
    """用户确认后的路由决策"""
    user_decision = state.get("user_decision", "").lower()
    
    if user_decision in ["e", "execute", "执行"]:
        print("✅ 开始执行分析")
        return "execute"
    elif user_decision in ["m", "modify", "修改"]:
        print("🔄 修改配置")
        return "modify"
    else:
        print("❌ 取消分析")
        return "cancel"