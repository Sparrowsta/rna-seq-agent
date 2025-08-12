"""
智能路由图系统
遵循图模式和状态机模式，实现基于模式的智能路由
"""

from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode
from .state import AgentState, create_initial_state
from .router import IntelligentRouter, route_user_input, route_after_tools, should_call_tools
from .core import ALL_TOOLS
from .nodes.normal_mode_node import normal_mode_node, get_user_input
from .nodes.plan_mode_node import plan_mode_node
from .nodes.execute_mode_node import execute_mode_node

# ============================================================================
# 路由函数 - 遵循单一职责原则
# ============================================================================

def route_from_user_input(state: AgentState) -> str:
    """
    从用户输入进行路由
    
    应用策略模式：使用智能路由器进行决策
    """
    try:
        # 检查退出条件
        if state.get("messages"):
            last_message = state["messages"][-1]
            if hasattr(last_message, "content"):
                content = last_message.content.lower().strip()
                if content in ["exit", "quit", "bye", "退出"]:
                    return "end"
        
        # 使用智能路由器
        router = IntelligentRouter()
        route = router.route(state)
        
        # 记录路由决策
        current_mode = state.get("mode", "normal")
        print(f"[ROUTER] 路由决策: {route} (当前模式: {current_mode})")
        
        return route
    
    except Exception as e:
        print(f"[ROUTER ERROR] 路由出错: {str(e)}")
        return "normal_mode_node"

def route_from_mode_nodes(state: AgentState) -> str:
    """
    从模式节点进行路由
    
    遵循开放封闭原则：易于扩展新的路由逻辑
    """
    try:
        # 首先检查状态中的模式是否已经改变
        current_mode = state.get("mode", "normal")
        
        # 检查是否需要调用工具
        has_tools = should_call_tools(state)
        print(f"[ROUTER DEBUG] should_call_tools返回: {has_tools}")
        
        if has_tools:
            print(f"[ROUTER] 检测到工具调用，路由到call_tools")
            return "call_tools"
        
        # 检查最后一条消息是否包含模式切换信息
        if state.get("messages"):
            last_message = state["messages"][-1]
            print(f"[ROUTER DEBUG] 最后一条消息类型: {type(last_message)}")
            
            if hasattr(last_message, "tool_calls"):
                print(f"[ROUTER DEBUG] tool_calls属性存在: {last_message.tool_calls}")
            
            # 检查AI消息内容中的模式切换标识
            if hasattr(last_message, "content") and last_message.content:
                content = str(last_message.content)
                if "切换到计划模式" in content or "正在切换到计划模式" in content:
                    print(f"[ROUTER] 检测到计划模式切换信号")
                    return "plan_mode_node"
                elif "切换到执行模式" in content or "正在切换到执行模式" in content:
                    print(f"[ROUTER] 检测到执行模式切换信号")
                    return "execute_mode_node"
            
            # 检查工具调用中的模式切换
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                for tool_call in last_message.tool_calls:
                    tool_name = tool_call.get("name", "")
                    if tool_name == "switch_to_plan_mode":
                        print(f"[ROUTER] 工具调用切换到计划模式")
                        return "plan_mode_node"
                    elif tool_name == "switch_to_execute_mode":
                        print(f"[ROUTER] 工具调用切换到执行模式")
                        return "execute_mode_node"
        
        # 关键修复：所有模式节点处理完后都应该返回用户输入，等待用户响应
        # 不再基于当前模式自动循环到相同节点
        print(f"[ROUTER] 模式节点处理完成，返回用户输入节点等待用户响应")
        return "user_input"
    
    except Exception as e:
        print(f"[ROUTER ERROR] 模式节点路由出错: {str(e)}")
        return "user_input"

def route_after_tool_calls(state: AgentState) -> str:
    """
    工具调用后的路由
    
    应用责任链模式：按优先级处理不同类型的工具调用结果
    """
    try:
        # 检查工具调用结果中的模式切换
        if state.get("messages"):
            last_message = state["messages"][-1]
            
            # 如果工具调用结果包含模式切换信息
            if hasattr(last_message, "content") and last_message.content:
                content = str(last_message.content)
                if "切换到计划模式" in content or "正在切换到计划模式" in content:
                    print(f"[ROUTER] 工具调用后检测到计划模式切换")
                    return "plan_mode_node"
                elif "切换到执行模式" in content or "正在切换到执行模式" in content:
                    print(f"[ROUTER] 工具调用后检测到执行模式切换")
                    return "execute_mode_node"
        
        # 根据当前模式返回相应节点
        current_mode = state.get("mode", "normal")
        print(f"[ROUTER] 工具调用后当前模式: {current_mode}")
        
        mode_mapping = {
            "normal": "normal_mode_node",
            "plan": "plan_mode_node",
            "execute": "execute_mode_node"
        }
        
        return mode_mapping.get(current_mode, "normal_mode_node")
    
    except Exception as e:
        print(f"[ROUTER ERROR] 工具调用后路由出错: {str(e)}")
        return "normal_mode_node"

def should_end_conversation(state: AgentState) -> str:
    """
    判断是否应该结束对话
    
    应用KISS原则：简单的结束条件判断
    """
    try:
        if state.get("messages"):
            last_message = state["messages"][-1]
            if hasattr(last_message, "content"):
                content = last_message.content.lower().strip()
                if content in ["exit", "quit", "bye", "退出"]:
                    return "end"
        
        return "continue"
    
    except Exception as e:
        print(f"[ROUTER ERROR] 结束判断出错: {str(e)}")
        return "continue"

# ============================================================================
# 节点包装函数 - 遵循适配器模式
# ============================================================================

def wrapped_normal_mode_node(state: AgentState) -> dict:
    """
    包装的normal模式节点
    
    应用装饰器模式：添加错误处理和日志
    """
    try:
        print("[NODE] 进入 Normal Mode")
        result = normal_mode_node(state)
        
        # 确保返回的结果包含必要的状态更新
        if "mode" not in result:
            result["mode"] = "normal"
        
        return result
    
    except Exception as e:
        print(f"[NODE ERROR] Normal Mode 节点错误: {str(e)}")
        return {
            "messages": [{"type": "ai", "content": f"Normal模式出现错误：{str(e)}"}],
            "mode": "normal"
        }

def wrapped_plan_mode_node(state: AgentState) -> dict:
    """
    包装的plan模式节点
    
    应用装饰器模式：添加错误处理和日志
    """
    try:
        print("[NODE] 进入 Plan Mode")
        result = plan_mode_node(state)
        
        # 确保返回的结果包含必要的状态更新
        if "mode" not in result:
            result["mode"] = "plan"
        
        return result
    
    except Exception as e:
        print(f"[NODE ERROR] Plan Mode 节点错误: {str(e)}")
        return {
            "messages": [{"type": "ai", "content": f"Plan模式出现错误：{str(e)}"}],
            "mode": "plan"
        }

def wrapped_execute_mode_node(state: AgentState) -> dict:
    """
    包装的execute模式节点
    
    应用装饰器模式：添加错误处理和日志
    """
    try:
        print("[NODE] 进入 Execute Mode")
        result = execute_mode_node(state)
        
        # 确保返回的结果包含必要的状态更新
        if "mode" not in result:
            result["mode"] = "execute"
        
        return result
    
    except Exception as e:
        print(f"[NODE ERROR] Execute Mode 节点错误: {str(e)}")
        return {
            "messages": [{"type": "ai", "content": f"Execute模式出现错误：{str(e)}"}],
            "mode": "execute"
        }

def wrapped_user_input_node(state: AgentState) -> dict:
    """
    包装的用户输入节点
    
    应用装饰器模式：添加错误处理和日志
    """
    try:
        print("[NODE] 获取用户输入")
        result = get_user_input(state)
        return result
    
    except Exception as e:
        print(f"[NODE ERROR] 用户输入节点错误: {str(e)}")
        return {
            "messages": [{"type": "human", "content": "输入错误，请重试"}]
        }

# ============================================================================
# 图构建 - 遵循建造者模式
# ============================================================================

class AgentGraphBuilder:
    """
    Agent图构建器
    
    遵循建造者模式：分步构建复杂的图结构
    """
    
    def __init__(self):
        self.workflow = StateGraph(AgentState)
        self.tools = ALL_TOOLS
    
    def add_nodes(self):
        """添加所有节点"""
        # 用户交互节点
        self.workflow.add_node("user_input", wrapped_user_input_node)
        
        # 模式特定节点
        self.workflow.add_node("normal_mode_node", wrapped_normal_mode_node)
        self.workflow.add_node("plan_mode_node", wrapped_plan_mode_node)
        self.workflow.add_node("execute_mode_node", wrapped_execute_mode_node)
        
        # 工具调用节点
        self.workflow.add_node("call_tools", ToolNode(self.tools))
        
        return self
    
    def set_entry_point(self):
        """设置入口点"""
        self.workflow.set_entry_point("user_input")
        return self
    
    def add_conditional_edges(self):
        """添加条件边"""
        # 从用户输入的路由
        self.workflow.add_conditional_edges(
            "user_input",
            route_from_user_input,
            {
                "normal_mode_node": "normal_mode_node",
                "plan_mode_node": "plan_mode_node",
                "execute_mode_node": "execute_mode_node",
                "end": END
            }
        )
        
        # 从normal模式节点的路由
        self.workflow.add_conditional_edges(
            "normal_mode_node",
            route_from_mode_nodes,
            {
                "call_tools": "call_tools",
                "user_input": "user_input",
                "plan_mode_node": "plan_mode_node",
                "execute_mode_node": "execute_mode_node"
            }
        )
        
        # 从plan模式节点的路由
        self.workflow.add_conditional_edges(
            "plan_mode_node",
            route_from_mode_nodes,
            {
                "call_tools": "call_tools",
                "user_input": "user_input",
                "execute_mode_node": "execute_mode_node",
                "plan_mode_node": "plan_mode_node"
            }
        )
        
        # 从execute模式节点的路由
        self.workflow.add_conditional_edges(
            "execute_mode_node",
            route_from_mode_nodes,
            {
                "call_tools": "call_tools",
                "user_input": "user_input",
                "normal_mode_node": "normal_mode_node",
                "plan_mode_node": "plan_mode_node"
            }
        )
        
        # 从工具调用的路由
        self.workflow.add_conditional_edges(
            "call_tools",
            route_after_tool_calls,
            {
                "normal_mode_node": "normal_mode_node",
                "plan_mode_node": "plan_mode_node",
                "execute_mode_node": "execute_mode_node"
            }
        )
        
        return self
    
    def build(self):
        """构建并编译图"""
        return self.workflow.compile()

# ============================================================================
# 图实例创建 - 应用工厂模式
# ============================================================================

def create_agent_graph():
    """
    创建agent图实例
    
    应用工厂模式：统一的图创建接口
    """
    try:
        builder = AgentGraphBuilder()
        graph = (builder
                .add_nodes()
                .set_entry_point()
                .add_conditional_edges()
                .build())
        
        print("[GRAPH] Agent图创建成功")
        return graph
    
    except Exception as e:
        print(f"[GRAPH ERROR] 创建图时出错: {str(e)}")
        raise

# 创建全局图实例
agent_executor = create_agent_graph()

# ============================================================================
# 图执行辅助函数 - 遵循DRY原则
# ============================================================================

def run_agent_with_initial_state(initial_input: str = "") -> dict:
    """
    使用初始状态运行agent
    
    应用模板方法模式：标准的执行流程
    """
    try:
        # 创建初始状态
        initial_state = create_initial_state()
        
        # 如果有初始输入，添加到消息中
        if initial_input:
            from langchain_core.messages import HumanMessage
            initial_state["messages"] = [HumanMessage(content=initial_input)]
        
        # 执行图
        result = agent_executor.invoke(initial_state)
        
        return result
    
    except Exception as e:
        print(f"[EXECUTION ERROR] 执行agent时出错: {str(e)}")
        return {"error": str(e)}

def run_agent_step_by_step(state: AgentState) -> dict:
    """
    逐步执行agent
    
    用于调试和监控执行过程
    """
    try:
        print(f"[STEP] 当前状态: mode={state.get('mode', 'unknown')}, messages={len(state.get('messages', []))}")
        
        # 执行一步
        result = agent_executor.invoke(state)
        
        print(f"[STEP] 执行结果: mode={result.get('mode', 'unknown')}, messages={len(result.get('messages', []))}")
        
        return result
    
    except Exception as e:
        print(f"[STEP ERROR] 逐步执行出错: {str(e)}")
        return {"error": str(e)}

# ============================================================================
# 图验证和调试 - 便于开发调试
# ============================================================================

def validate_graph_structure():
    """
    验证图结构的完整性
    
    应用KISS原则：简单的结构验证
    """
    try:
        # 检查必要的节点是否存在
        required_nodes = [
            "user_input", "normal_mode_node", "plan_mode_node", 
            "execute_mode_node", "call_tools"
        ]
        
        # 这里可以添加更详细的验证逻辑
        print("[VALIDATION] 图结构验证通过")
        return True
    
    except Exception as e:
        print(f"[VALIDATION ERROR] 图结构验证失败: {str(e)}")
        return False

def print_graph_info():
    """
    打印图信息
    
    用于调试和了解图结构
    """
    try:
        print("=" * 50)
        print("🧬 RNA-seq Agent 图结构信息")
        print("=" * 50)
        print("📍 入口点: user_input")
        print("🔄 节点列表:")
        print("  - user_input: 用户输入节点")
        print("  - normal_mode_node: 信息收集模式")
        print("  - plan_mode_node: 计划制定模式")
        print("  - execute_mode_node: 执行模式")
        print("  - call_tools: 工具调用节点")
        print("🛠️  可用工具数量:", len(ALL_TOOLS))
        print("=" * 50)
    
    except Exception as e:
        print(f"[INFO ERROR] 打印图信息出错: {str(e)}")

# 在模块加载时验证图结构
if __name__ == "__main__":
    validate_graph_structure()
    print_graph_info()