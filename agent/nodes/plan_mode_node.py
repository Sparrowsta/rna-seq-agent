"""
Plan Mode节点 - 制定分析计划和修改nextflow参数
遵循单一职责原则：专门处理plan模式下的计划制定和参数配置
采用JSON-first架构，与normal模式保持一致
简化版，移除execution_phase复杂度
"""

import logging
from typing import Dict, Any, List, Tuple
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_state_mode
from ..core import create_structured_chain_for_mode
from ..ui_manager import get_ui_manager

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PlanModeHandler:
    """
    Plan模式处理器
    
    遵循单一职责原则：专门处理plan模式的业务逻辑
    采用JSON-first架构，与normal模式保持一致
    """
    
    def __init__(self):
        # 使用Pydantic结构化输出链
        self.chain = create_structured_chain_for_mode("plan")
    
    def _process_llm_response(self, response) -> Tuple[AIMessage, List[Dict[str, Any]]]:
        """
        处理LLM的结构化响应（.with_structured_output()返回Pydantic模型实例）
        DeepSeek + json_mode确保总是返回结构化数据，无需手动JSON解析
        
        返回: (AIMessage, tool_calls列表)
        """
        try:
            # 调试日志：查看响应类型和内容
            logger.info(f"Plan模式收到响应类型: {type(response)}")
            logger.info(f"Plan模式收到响应内容: {repr(response)[:200]}...")
            
            # 处理情况1：Pydantic模型实例（期望的格式）
            if hasattr(response, 'reasoning') and hasattr(response, 'plan_steps'):  # PlanModeResponse模型
                logger.info(f"Plan模式收到Pydantic模型响应")
                
                # 提取响应信息
                reasoning = getattr(response, "reasoning", "计划制定完成")
                next_action = getattr(response, "next_action", "")
                plan_steps = getattr(response, "plan_steps", [])
                config_changes = getattr(response, "config_changes", {})
                ready_to_execute = getattr(response, "ready_to_execute", False)
                tool_calls = getattr(response, "tool_calls", [])
                
                # 构建详细响应 - 使用reasoning作为主要回复内容
                detailed_response = reasoning
                if next_action:
                    detailed_response += f"\n\n➡️ **下一步行动**: {next_action}"
                if plan_steps:
                    detailed_response += "\n\n📋 **分析计划步骤：**\n"
                    detailed_response += "\n".join([f"  {step}" for step in plan_steps])
                
                if config_changes:
                    detailed_response += "\n\n⚙️ **配置更新：**\n"
                    for key, value in config_changes.items():
                        detailed_response += f"  - {key}: {value}\n"
                
                # 如果已准备好执行，添加执行提示
                if ready_to_execute:
                    detailed_response += "\n\n🚀 **配置完成！**\n"
                    detailed_response += "所有参数已配置完成，可以开始执行RNA-seq分析。\n"
                    detailed_response += "请输入 `/execute` 或 `/开始执行` 开始分析流程。"
                
                logger.info(f"Plan模式提取到 {len(tool_calls)} 个工具调用")
                
                # 创建AIMessage
                ai_message = AIMessage(content=detailed_response)
                
                # 如果有工具调用，设置为消息的tool_calls属性
                if tool_calls:
                    langchain_tool_calls = []
                    for i, tool_call in enumerate(tool_calls):
                        # 处理Pydantic模型中的工具调用
                        if hasattr(tool_call, 'tool_name'):
                            # tool_call是ToolCall Pydantic模型实例
                            tool_call_obj = {
                                "name": tool_call.tool_name,
                                "args": tool_call.parameters,
                                "id": f"call_plan_{i}",
                                "type": "tool_call"
                            }
                        else:
                            # tool_call是字典格式
                            tool_call_obj = {
                                "name": tool_call.get("tool_name"),
                                "args": tool_call.get("parameters", {}),
                                "id": f"call_plan_{i}",
                                "type": "tool_call"
                            }
                        langchain_tool_calls.append(tool_call_obj)
                    
                    ai_message.tool_calls = langchain_tool_calls
                
                return ai_message, tool_calls
            
            # 处理情况2：纯字符串回复
            elif isinstance(response, str):
                logger.info("Plan模式收到字符串响应")
                ai_message = AIMessage(content=response)
                return ai_message, []
            
            else:
                # 未知格式，尝试字符串化
                logger.warning(f"Plan模式收到未知格式响应: {type(response)}")
                ai_message = AIMessage(content=str(response))
                return ai_message, []
                
        except Exception as e:
            logger.error(f"Plan模式处理响应时出错: {str(e)}")
            error_message = AIMessage(content=f"处理响应时出错: {str(e)}")
            return error_message, []
    
    def handle_plan_modification(self, state: AgentState) -> Dict[str, Any]:
        """
        处理计划修改请求
        
        应用单一职责原则：专门处理用户的计划修改需求
        """
        try:
            # 检查是否有最新的用户输入
            messages = state.get("messages", [])
            if not messages:
                return {"messages": [HumanMessage(content="未检测到用户输入，请提供分析需求或计划修改要求")]}
            
            logger.info("Plan模式开始处理计划修改...")
            
            # 调用LLM进行计划分析和制定
            response = self.chain.invoke(state)
            ai_message, tool_calls = self._process_llm_response(response)
            
            # 返回结果
            result = {"messages": [ai_message]}
            
            # 如果有配置更新，同时更新nextflow_config
            if hasattr(response, 'config_changes') and response.config_changes:
                logger.info(f"更新nextflow配置: {response.config_changes}")
                # 更新配置
                updated_config = state.get("nextflow_config", {}).copy()
                updated_config.update(response.config_changes)
                result["nextflow_config"] = updated_config
            
            return result
            
        except Exception as e:
            logger.error(f"Plan模式处理计划修改时出错: {str(e)}")
            error_message = HumanMessage(content=f"处理计划修改时出错: {str(e)}")
            return {"messages": [error_message]}
    
    def handle_mode_switch_request(self, state: AgentState) -> Dict[str, Any]:
        """
        检查并处理模式切换请求
        
        优先级最高的处理逻辑，检查工具调用中的模式切换指令
        """
        try:
            messages = state.get("messages", [])
            if not messages:
                return None
            
            # 检查最后一条消息是否包含工具调用
            last_message = messages[-1]
            
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                for tool_call in last_message.tool_calls:
                    tool_name = tool_call.get("name", "")
                    
                    # 检查执行模式切换
                    if tool_name == "switch_to_execute_mode":
                        logger.info("检测到切换到执行模式的工具调用")
                        return {"mode": "execute"}
                    
                    # 检查普通模式切换
                    elif tool_name == "switch_to_normal_mode":
                        logger.info("检测到切换到普通模式的工具调用")
                        return {"mode": "normal"}
            
            return None
            
        except Exception as e:
            logger.error(f"Plan模式处理模式切换时出错: {str(e)}")
            return None
    
    def _auto_initialize_agent_state(self, state: AgentState) -> Dict[str, Any]:
        """
        自动初始化AgentState
        
        为首次进入Plan模式的用户自动收集基本信息并制定初始计划
        """
        try:
            logger.info("开始自动初始化AgentState...")
            
            # 创建初始化消息，指导AI收集信息并制定计划
            init_message = HumanMessage(
                content="请帮我制定RNA-seq分析计划。请先查看当前环境中的FASTQ文件和可用的基因组配置，然后制定详细的分析计划。"
            )
            
            # 将初始化消息添加到state中
            updated_state = state.copy()
            updated_messages = state.get("messages", []) + [init_message]
            updated_state["messages"] = updated_messages
            
            # 调用LLM进行初始化
            response = self.chain.invoke(updated_state)
            ai_message, tool_calls = self._process_llm_response(response)
            
            # 构建初始化结果
            result = {"messages": updated_messages + [ai_message]}
            
            # 如果有配置更新，更新nextflow_config
            if hasattr(response, 'config_changes') and response.config_changes:
                logger.info(f"初始化时更新nextflow配置: {response.config_changes}")
                updated_config = state.get("nextflow_config", {}).copy()
                updated_config.update(response.config_changes)
                result["nextflow_config"] = updated_config
            
            return result
            
        except Exception as e:
            logger.error(f"自动初始化AgentState时出错: {str(e)}")
            return None

def plan_mode_node(state: AgentState) -> Dict[str, Any]:
    """
    Plan模式主节点函数 - 简化版
    
    直接处理用户的计划制定和修改请求
    """
    logger.info("进入Plan模式节点")
    
    try:
        # 获取UI管理器
        ui_manager = get_ui_manager()
        current_mode = state.get("mode", "normal")
        
        logger.info(f"Plan模式当前模式: {current_mode}")
        
        # 创建处理器
        handler = PlanModeHandler()
        
        # 检查是否需要处理模式切换（优先级最高）
        mode_switch_result = handler.handle_mode_switch_request(state)
        if mode_switch_result:
            logger.info("检测到模式切换请求")
            if mode_switch_result.get("mode") == "execute":
                ui_manager.show_mode_switch("plan", "execute", "计划已确认，开始执行")
            return mode_switch_result
        
        # 首次进入Plan模式 - 初始化
        if current_mode != "plan":
            logger.info("首次进入Plan模式，进行初始化...")
            ui_manager.show_mode_switch(current_mode, "plan", "开始制定分析计划")
            
            # 自动初始化AgentState
            plan_response = handler._auto_initialize_agent_state(state)
            if plan_response:
                logger.info("AgentState初始化完成")
                plan_response["mode"] = "plan"
                return plan_response
            
            # 如果自动初始化失败，返回错误信息
            return {
                "messages": [HumanMessage(content="自动初始化失败，请手动提供分析需求")],
                "mode": "plan"
            }
        
        # Plan模式下的用户交互处理
        result = handler.handle_plan_modification(state)
        result["mode"] = "plan"
        return result
        
    except Exception as e:
        logger.error(f"Plan模式节点执行出错: {str(e)}")
        error_message = HumanMessage(content=f"执行过程中发生错误: {str(e)}")
        return {
            "messages": [error_message],
            "mode": "plan"
        }