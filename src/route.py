from langgraph.graph import END
from .state import PlanExecuteState

async def should_continue(state: PlanExecuteState) -> str:
    """决定是否继续执行"""
    plan = state.get("plan", [])
    if plan:
        return "execute"
    else:
        return END