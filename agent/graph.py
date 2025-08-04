from langgraph.graph import END, StateGraph
from .state import AgentState
from .nodes.chat_nodes import call_model, get_user_input

def should_continue(state: AgentState):
    last_message = state["messages"][-1]
    if last_message.content.lower() == "exit":
        return "end"
    return "continue"

workflow = StateGraph(AgentState)

workflow.add_node("user_input", get_user_input)
workflow.add_node("agent", call_model)

workflow.set_entry_point("user_input")

workflow.add_conditional_edges(
    "user_input",
    should_continue,
    {
        "continue": "agent",
        "end": END,
    },
)
workflow.add_edge("agent", "user_input")

agent_executor = workflow.compile()