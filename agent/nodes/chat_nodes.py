from langchain_core.messages import HumanMessage
from ..state import AgentState
from ..core import prompt, llm_with_tools


def get_user_input(state: AgentState):
    last_message = state["messages"][-1] if state["messages"] else None
    if last_message and last_message.type != "human":
        print(f"AI: {last_message.content}")

    user_input = input("You: ")
    return {"messages": [HumanMessage(content=user_input)]}


def call_model(state: AgentState):
    chain = prompt | llm_with_tools
    response = chain.invoke(
        {"messages": state["messages"]},
        {"recursion_limit": 50}
    )
    return {"messages": [response]}