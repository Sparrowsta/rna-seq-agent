from langgraph.graph import END
from .state import NormalNodeState, PrepareNodeState, ReplanNodeState, UserConfirmState

async def route_from_normal(state: NormalNodeState) -> str:
    """Normal节点后的路由决策"""
    routing_decision = state.get("routing_decision", "normal")
    
    if routing_decision == "plan":
        print("➡️ 用户提出分析需求，进入Plan流程")
        return "plan"
    else:
        print("🔄 继续Normal模式交互")
        return "normal"

async def should_continue(state: PrepareNodeState) -> str:
    """决定是否继续执行"""
    plan = state.get("plan", [])
    if plan:
        return "execute"
    else:
        return END

async def route_after_confirm(state: UserConfirmState) -> str:
    """用户确认后的路由决策"""
    user_decision = state.get("user_decision", "").lower()
    
    if user_decision in ["e", "execute", "执行"]:
        print("✅ 用户选择执行分析，结束配置流程")
        return "execute"
    elif user_decision in ["m", "modify", "修改"]:
        print("🔄 用户选择修改配置，进入Replan流程")
        return "modify"
    else:  # cancel or 取消 or empty
        print("❌ 用户取消分析或输入为空，结束流程")
        return "cancel"