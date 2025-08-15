"""
智能路由系统
遵循策略模式和单一职责原则，处理基于状态和工具调用的路由决策
"""

from typing import Dict, Any, Optional
from .state import AgentState

class RouterStrategy:
    """
    路由策略基类
    
    遵循开放封闭原则：可扩展新的路由策略而不修改现有代码
    """
    
    def route(self, state: AgentState) -> str:
        """路由决策方法，子类必须实现"""
        raise NotImplementedError

class ToolCallRouter(RouterStrategy):
    """
    基于工具调用的路由策略
    
    应用KISS原则：简单直接的工具调用检测
    """
    
    def __init__(self):
        # 工具调用到目标节点的映射
        self.tool_to_node_mapping = {
            "switch_to_plan_mode": "plan_mode_node",
            "switch_to_execute_mode": "execute_mode_node",
            "execute_nextflow_pipeline": "execute_mode_node",
        }
    
    def route(self, state: AgentState) -> Optional[str]:
        """基于工具调用进行路由"""
        if not state.get("messages"):
            return None
        
        last_message = state["messages"][-1]
        
        # 检查是否有工具调用
        if hasattr(last_message, "tool_calls") and last_message.tool_calls:
            # 查找模式切换相关的工具调用
            for tool_call in last_message.tool_calls:
                tool_name = tool_call.get("name", "")
                if tool_name in self.tool_to_node_mapping:
                    return self.tool_to_node_mapping[tool_name]
        
        # 检查工具调用结果中的模式切换信息
        if hasattr(last_message, "content") and last_message.content:
            content = str(last_message.content)
            if "切换到计划模式" in content or "正在切换到计划模式" in content:
                return "plan_mode_node"
            elif "切换到执行模式" in content or "正在切换到执行模式" in content:
                return "execute_mode_node"
        
        return None

class ModeBasedRouter(RouterStrategy):
    """
    基于当前模式的路由策略
    
    遵循单一职责原则：专门处理模式到节点的映射
    """
    
    def __init__(self):
        # 模式到节点的映射
        self.mode_to_node_mapping = {
            "normal": "normal_mode_node",
            "plan": "plan_mode_node",
            "execute": "execute_mode_node"
        }
    
    def route(self, state: AgentState) -> str:
        """基于当前模式进行路由"""
        current_mode = state.get("mode", "normal")
        return self.mode_to_node_mapping.get(current_mode, "normal_mode_node")

class ConversationRouter(RouterStrategy):
    """
    基于对话内容的路由策略
    
    应用YAGNI原则：只实现当前需要的基础对话路由
    """
    
    def route(self, state: AgentState) -> Optional[str]:
        """基于对话内容进行路由"""
        if not state.get("messages"):
            return None
        
        last_message = state["messages"][-1]
        
        # 检查是否为用户输入
        if hasattr(last_message, "type") and last_message.type == "human":
            content = last_message.content.lower().strip()
            
            # 检查退出指令
            if content in ["exit", "quit", "bye", "退出"]:
                return "end"
            
            # 关键修复：检查特殊命令（最高优先级）
            if content in ["/plan", "/开始计划", "/制定计划"]:
                print(f"[CONVERSATION ROUTER] 检测到计划模式命令: {content}")
                return "plan_mode_node"
            elif content in ["/execute", "/开始执行", "/执行"]:
                print(f"[CONVERSATION ROUTER] 检测到执行模式命令: {content}")
                return "execute_mode_node"
            
            # 检查是否需要工具调用
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                return "call_tools"
        
        return None

class IntelligentRouter:
    """
    智能路由器主类
    
    遵循组合模式：组合多种路由策略
    应用责任链模式：按优先级尝试不同的路由策略
    """
    
    def __init__(self):
        # 按优先级排序的路由策略
        self.strategies = [
            ConversationRouter(),  # 最高优先级：对话控制
            ToolCallRouter(),      # 中等优先级：工具调用
            ModeBasedRouter()      # 最低优先级：默认模式路由
        ]
    
    def route(self, state: AgentState) -> str:
        """
        智能路由决策
        
        应用责任链模式：依次尝试各种路由策略
        """
        # 按优先级尝试各种路由策略
        for strategy in self.strategies:
            result = strategy.route(state)
            if result:
                return result
        
        # 默认路由到normal模式
        return "normal_mode_node"

# ============================================================================
# 路由辅助函数 - 遵循DRY原则
# ============================================================================

def route_after_tools(state: AgentState) -> str:
    """
    工具调用后的路由决策
    
    直接根据模式进行路由，无需复杂的执行阶段判断
    """
    current_mode = state.get("mode", "normal")
    
    print(f"[ROUTER] 工具调用后 - 模式: {current_mode}")
    
    # 直接根据模式进行路由
    mode_mapping = {
        "normal": "normal_mode_node", 
        "plan": "plan_mode_node",
        "execute": "execute_mode_node"
    }
    
    next_node = mode_mapping.get(current_mode, "normal_mode_node")
    print(f"[ROUTER] 路由到: {next_node}")
    return next_node

def route_user_input(state: AgentState) -> str:
    """
    用户输入后的路由决策
    
    遵循单一职责原则：专门处理用户输入路由
    """
    if not state.get("messages"):
        return "normal_mode_node"
    
    last_message = state["messages"][-1]
    
    # 检查退出指令
    if hasattr(last_message, "content"):
        content = last_message.content.lower().strip()
        if content in ["exit", "quit", "bye", "退出"]:
            return "end"
    
    # 使用智能路由器进行决策
    router = IntelligentRouter()
    return router.route(state)

def should_call_tools(state: AgentState) -> bool:
    """
    判断是否需要调用工具
    
    所有模式都支持工具调用
    """
    if not state.get("messages"):
        print("[ROUTER DEBUG] should_call_tools: 没有消息")
        return False
    
    last_message = state["messages"][-1]
    print(f"[ROUTER DEBUG] should_call_tools: 最后一条消息类型: {type(last_message)}")
    print(f"[ROUTER DEBUG] should_call_tools: 消息有tool_calls属性吗? {hasattr(last_message, 'tool_calls')}")
    
    if hasattr(last_message, "tool_calls"):
        tool_calls = last_message.tool_calls
        print(f"[ROUTER DEBUG] should_call_tools: tool_calls值: {tool_calls}")
        print(f"[ROUTER DEBUG] should_call_tools: tool_calls长度: {len(tool_calls) if tool_calls else 0}")
        
        result = (tool_calls and len(tool_calls) > 0)
        print(f"[ROUTER DEBUG] should_call_tools: 最终结果: {result}")
        return result
    else:
        print("[ROUTER DEBUG] should_call_tools: 消息没有tool_calls属性")
        return False

def get_next_node_after_mode_switch(state: AgentState, switched_mode: str) -> str:
    """
    模式切换后的下一个节点决策
    
    遵循开放封闭原则：易于扩展新的模式处理
    """
    mode_mapping = {
        "plan": "plan_mode_node",
        "execute": "execute_mode_node",
        "normal": "normal_mode_node"
    }
    
    return mode_mapping.get(switched_mode, "normal_mode_node")

# ============================================================================
# 路由配置和验证 - 遵循配置分离原则
# ============================================================================

class RouteConfig:
    """
    路由配置类
    
    遵循单一职责原则：专门管理路由配置
    """
    
    # 有效的节点名称
    VALID_NODES = {
        "normal_mode_node",
        "plan_mode_node", 
        "execute_mode_node",
        "call_tools",
        "user_input",
        "end"
    }
    
    # 有效的模式
    VALID_MODES = {
        "normal",
        "plan", 
        "execute"
    }
    
    @classmethod
    def validate_route(cls, route: str) -> bool:
        """验证路由是否有效"""
        return route in cls.VALID_NODES
    
    @classmethod
    def validate_mode(cls, mode: str) -> bool:
        """验证模式是否有效"""
        return mode in cls.VALID_MODES

def create_router() -> IntelligentRouter:
    """
    创建路由器实例
    
    应用工厂模式：统一的路由器创建
    """
    return IntelligentRouter()

# ============================================================================
# 路由调试和日志 - 便于开发调试
# ============================================================================

def log_routing_decision(state: AgentState, route: str, reason: str = "") -> None:
    """
    记录路由决策日志
    
    应用YAGNI原则：简单的日志记录，便于调试
    """
    current_mode = state.get("mode", "unknown")
    message_count = len(state.get("messages", []))
    
    log_message = f"路由决策: {route} (模式: {current_mode}, 消息数: {message_count})"
    if reason:
        log_message += f" - 原因: {reason}"
    
    # 这里可以集成到正式的日志系统
    print(f"[ROUTER] {log_message}")