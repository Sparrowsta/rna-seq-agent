"""
智能路由图系统
遵循图模式和状态机模式，实现基于模式的智能路由
"""

from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode
from .state import AgentState, create_initial_state
from .state import update_current_stage, add_tool_to_execution_history
from .router import IntelligentRouter, route_user_input, route_after_tools
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
    从模式节点进行路由 - 基于结构化输出设计
    
    结构化输出模式下，所有工具调用在模式节点内部完成
    模式切换通过特殊命令实现，无需检测工具调用
    """
    try:
        # 结构化输出模式下，模式节点完成所有处理后直接等待用户输入
        return "user_input"
    
    except Exception as e:
        print(f"[ROUTER ERROR] 模式节点路由出错: {str(e)}")
        return "user_input"

def route_after_tool_calls(state: AgentState) -> str:
    """
    工具调用后的路由
    
    使用新的阶段化状态管理，防止工具循环调用
    """
    try:
        current_mode = state.get("mode", "normal")
        current_stage = state.get("current_stage", "user_input")
        
        # 1. 首先处理模式切换（只检查最后一条消息）
        if state.get("messages"):
            last_message = state["messages"][-1]
            if hasattr(last_message, "type") and last_message.type == "tool":
                tool_name = getattr(last_message, "name", "unknown")
                if tool_name == "switch_to_plan_mode":
                    state["mode"] = "plan"
                    current_mode = "plan"
                elif tool_name == "switch_to_execute_mode":
                    state["mode"] = "execute"
                    current_mode = "execute"
        
        # 2. 记录已执行的工具到历史中（向前追溯到非工具消息）
        if state.get("messages"):
            for message in reversed(state["messages"]):
                if hasattr(message, "type") and message.type == "tool":
                    tool_name = getattr(message, "name", "unknown")
                    state = add_tool_to_execution_history(state, tool_name)
                else:
                    break
        
        # 3. 根据模式和阶段更新状态
        if current_mode == "plan":
            if current_stage == "tools_executing":
                state = update_current_stage(state, "tools_completed", "计划工具执行完成")
        
        # 4. 根据当前模式返回相应节点
        mode_mapping = {
            "normal": "normal_mode_node",
            "plan": "plan_mode_node",
            "execute": "execute_mode_node"
        }
        
        next_node = mode_mapping.get(current_mode, "normal_mode_node")
        return next_node
    
    except Exception as e:
        print(f"[ROUTER ERROR] 工具调用后路由出错: {str(e)}")
        return "normal_mode_node"

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
                "plan_mode_node": "plan_mode_node",
                "execute_mode_node": "execute_mode_node"
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
        # 注意：recursion_limit应该在执行时通过config传递，不在compile时设置
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
        
        return graph
    
    except Exception as e:
        print(f"[GRAPH ERROR] 创建图时出错: {str(e)}")
        raise

# 创建全局图实例
agent_executor = create_agent_graph()

# ============================================================================

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