from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode
from .state import AgentState
from .nodes.chat_nodes import call_model, get_user_input
from .core import tools

def route_user_input(state: AgentState):
    """Routes user input to agent or ends the conversation."""
    last_message = state["messages"][-1]
    if last_message.content.lower() == "exit":
        return "end"
    return "agent"

def route_after_agent(state: AgentState):
    """Checks for tool calls and routes accordingly."""
    last_message = state["messages"][-1]
    if hasattr(last_message, "tool_calls") and last_message.tool_calls:
        return "call_tools"
    return "continue"

workflow = StateGraph(AgentState)

workflow.add_node("user_input", get_user_input)
workflow.add_node("agent", call_model)
workflow.add_node("tools", ToolNode(tools))

workflow.set_entry_point("user_input")

workflow.add_conditional_edges(
    "user_input",
    route_user_input,
    {
        "agent": "agent",
        "end": END,
    },
)

workflow.add_conditional_edges(
    "agent",
    route_after_agent,
    {
        "continue": "user_input",
        "call_tools": "tools",
    },
)
workflow.add_edge("tools", "agent")

agent_executor = workflow.compile()