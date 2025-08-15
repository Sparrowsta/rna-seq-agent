"""
Normal Mode节点 - 信息收集和用户交互
遵循单一职责原则：专门处理normal模式下的用户交互和信息收集
"""

import logging
from typing import Dict, Any, List, Tuple
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_state_mode
from ..core import create_chain_for_mode, create_structured_chain_for_mode

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NormalModeHandler:
    """
    Normal模式处理器
    
    遵循单一职责原则：专门处理normal模式的业务逻辑
    """
    
    def __init__(self):
        # 使用不带工具的LLM链用于JSON格式输出
        self.chain = create_structured_chain_for_mode("normal")  
        self.structured_chain = create_structured_chain_for_mode("normal")
    
    def _process_llm_response(self, response) -> Tuple[AIMessage, List[Dict[str, Any]]]:
        """
        处理LLM的结构化响应（.with_structured_output()返回Pydantic模型实例）
        
        返回: (AIMessage, tool_calls列表)
        """
        try:
            # 调试日志：查看响应类型和内容
            logger.info(f"Normal模式收到响应类型: {type(response)}")
            logger.info(f"Normal模式收到响应内容: {response}")
            
            # .with_structured_output()返回Pydantic模型实例
            if hasattr(response, 'response') and hasattr(response, 'tool_calls'):  # NormalModeResponse模型
                logger.info(f"Normal模式收到Pydantic模型响应")
                
                # 提取响应信息
                user_message = getattr(response, "response", "响应完成")
                suggested_actions = getattr(response, "suggested_actions", [])
                need_more_info = getattr(response, "need_more_info", False)
                tool_calls = getattr(response, "tool_calls", [])
                
                # 构建详细响应
                detailed_response = user_message
                if suggested_actions:
                    detailed_response += "\n\n💡 **建议操作：**\n"
                    detailed_response += "\n".join([f"  - {action}" for action in suggested_actions])
                
                if need_more_info:
                    detailed_response += "\n\n❓ 需要更多信息才能继续。"
                
                logger.info(f"Normal模式提取到 {len(tool_calls)} 个工具调用")
                
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
                                "id": f"call_{i}",
                                "type": "tool_call"
                            }
                        else:
                            # tool_call是字典格式
                            tool_call_obj = {
                                "name": tool_call.get("tool_name"),
                                "args": tool_call.get("parameters", {}),
                                "id": f"call_{i}",
                                "type": "tool_call"
                            }
                        langchain_tool_calls.append(tool_call_obj)
                    
                    ai_message.tool_calls = langchain_tool_calls
                    logger.info(f"Normal模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            elif isinstance(response, str) and 'query_fastq_files' in response:
                # 特殊处理：如果DeepSeek直接返回了工具名字符串
                logger.warning(f"DeepSeek返回了工具名字符串而不是结构化响应: {response}")
                
                # 尝试手动构造工具调用
                ai_message = AIMessage(content="正在查询FASTQ文件信息...")
                tool_call_obj = {
                    "name": "query_fastq_files",
                    "args": {},
                    "id": "call_manual_0",
                    "type": "tool_call"
                }
                ai_message.tool_calls = [tool_call_obj]
                
                manual_tool_calls = [{
                    "tool_name": "query_fastq_files",
                    "parameters": {},
                    "reason": "用户请求查看FASTQ文件"
                }]
                
                logger.info(f"手动构造了工具调用: query_fastq_files")
                return ai_message, manual_tool_calls
            elif isinstance(response, dict):
                # 兼容旧的dict格式返回
                logger.info(f"Normal模式收到dict格式响应: {list(response.keys())}")
                
                # 提取响应信息
                user_message = response.get("response", "响应完成")
                suggested_actions = response.get("suggested_actions", [])
                need_more_info = response.get("need_more_info", False)
                tool_calls = response.get("tool_calls", [])
                
                # 构建详细响应
                detailed_response = user_message
                if suggested_actions:
                    detailed_response += "\n\n💡 **建议操作：**\n"
                    detailed_response += "\n".join([f"  - {action}" for action in suggested_actions])
                
                if need_more_info:
                    detailed_response += "\n\n❓ 需要更多信息才能继续。"
                
                logger.info(f"Normal模式提取到 {len(tool_calls)} 个工具调用")
                
                # 创建AIMessage
                ai_message = AIMessage(content=detailed_response)
                
                # 如果有工具调用，设置为消息的tool_calls属性
                if tool_calls:
                    langchain_tool_calls = []
                    for i, tool_call in enumerate(tool_calls):
                        tool_call_obj = {
                            "name": tool_call.get("tool_name"),
                            "args": tool_call.get("parameters", {}),
                            "id": f"call_{i}",
                            "type": "tool_call"
                        }
                        langchain_tool_calls.append(tool_call_obj)
                    
                    ai_message.tool_calls = langchain_tool_calls
                    logger.info(f"Normal模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            else:
                # 降级处理：如果不是期望的格式
                logger.warning(f"Normal模式收到未知响应格式: {type(response)}")
                content = str(response) if response else "响应为空"
                return AIMessage(content=content), []
            
        except Exception as e:
            logger.error(f"Normal模式处理响应时出错: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            return AIMessage(content=f"处理响应时出现错误: {str(e)}"), []
    
    def process_user_input(self, state: AgentState) -> Dict[str, Any]:
        """
        处理用户输入
        
        应用KISS原则：简单直接的用户输入处理
        """
        try:
            # 确保当前处于normal模式
            if state.get("mode") != "normal":
                logger.warning(f"Expected normal mode, but got {state.get('mode')}")
                state = update_state_mode(state, "normal")
            
            # 获取最后一条用户消息
            user_input = ""
            if state.get("messages"):
                last_message = state["messages"][-1]
                if hasattr(last_message, "content"):
                    user_input = last_message.content
            
            # 调用LLM处理用户输入，传入input参数
            try:
                logger.info(f"调用LLM链处理用户输入: {user_input}")
                response = self.chain.invoke({
                    "messages": state["messages"],
                    "input": user_input
                })
                logger.info(f"LLM链调用成功，响应类型: {type(response)}")
                logger.info(f"LLM响应内容: {repr(response)}")
            except Exception as chain_error:
                logger.error(f"LLM链调用失败: {str(chain_error)}")
                import traceback
                logger.error(f"链调用错误堆栈: {traceback.format_exc()}")
                raise chain_error
            
            # 处理LLM的结构化响应
            try:
                parsed_response, tool_calls = self._process_llm_response(response)
                logger.info(f"响应解析成功")
            except Exception as parse_error:
                logger.error(f"响应解析失败: {str(parse_error)}")
                import traceback
                logger.error(f"解析错误堆栈: {traceback.format_exc()}")
                raise parse_error
            
            # 如果有工具调用，添加到响应中
            if tool_calls:
                logger.info(f"检测到工具调用，准备返回带有tool_calls的消息")
                logger.info(f"parsed_response.tool_calls: {getattr(parsed_response, 'tool_calls', None)}")
            else:
                logger.info("没有检测到工具调用")
            
            logger.info(f"Normal mode response generated: {type(response)}")
            logger.info(f"返回的消息tool_calls属性: {hasattr(parsed_response, 'tool_calls')} - {getattr(parsed_response, 'tool_calls', None)}")
            
            return {"messages": [parsed_response]}
        
        except Exception as e:
            logger.error(f"Error in normal mode processing: {str(e)}")
            error_message = AIMessage(
                content=f"抱歉，处理您的请求时出现错误：{str(e)}。请重试或联系技术支持。"
            )
            return {"messages": [error_message]}
    
    def handle_mode_switch_request(self, state: AgentState) -> Dict[str, Any]:
        """
        处理模式切换请求
        
        遵循开放封闭原则：易于扩展新的模式切换逻辑
        只通过工具调用进行模式切换，不再使用关键词检测
        """
        try:
            # 检查最后一条消息是否包含模式切换工具调用
            if not state.get("messages"):
                return {}
                
            last_message = state["messages"][-1]
            
            # 只检查工具调用，移除关键词检测
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                for tool_call in last_message.tool_calls:
                    if tool_call.get("name") == "switch_to_plan_mode":
                        logger.info("Switching to plan mode requested via tool call")
                        # 保持完整状态，只更新模式
                        result = dict(state)  # 复制现有状态
                        result["mode"] = "plan"
                        result["messages"] = state["messages"] + [
                            AIMessage(content="🔄 已切换到计划模式！现在开始制定RNA-seq分析计划...")
                        ]
                        return result
            
            # 移除基于消息内容的模式切换检测
            # 所有模式切换都应通过LLM的工具调用来实现
            
            return {}
        
        except Exception as e:
            logger.error(f"Error handling mode switch: {str(e)}")
            return {}
    
    def provide_guidance(self, state: AgentState) -> Dict[str, Any]:
        """
        提供用户指导
        
        应用YAGNI原则：只提供当前需要的基础指导
        """
        guidance_message = AIMessage(content="""
欢迎使用RNA-seq智能助手！我可以帮助您：

🔍 **信息查询**：
- 查看目录内容和FASTQ文件
- 查询可用的基因组配置
- 了解当前的处理配置

📋 **处理准备**：
- 回答RNA-seq相关问题
- 帮助您准备处理所需的文件
- 制定个性化的处理计划

🚀 **开始处理**：
当您准备好开始处理时，只需告诉我"开始分析"或"制定计划"

请告诉我您想了解什么，或者您有什么FASTQ文件需要处理？
        """)
        
        return {"messages": [guidance_message]}

def get_user_input(state: AgentState) -> Dict[str, Any]:
    """
    获取用户输入的节点函数 - 使用增强的UI管理器
    
    遵循单一职责原则：专门处理用户输入获取
    """
    try:
        from ..ui_manager import get_ui_manager
        
        ui_manager = get_ui_manager()
        
        # 显示AI的最后回复（如果有）
        if state.get("messages"):
            last_message = state["messages"][-1]
            if hasattr(last_message, "type") and last_message.type == "ai":
                ui_manager.show_ai_response(last_message.content, "normal")
        
        # 获取用户输入
        user_input = ui_manager.get_user_input("请告诉我您的需求", "normal")
        
        if not user_input or user_input.lower() in ["exit", "quit", "退出"]:
            return {"messages": [HumanMessage(content="exit")]}
        
        logger.info(f"User input received: {user_input[:50]}...")
        
        # 更新状态以包含用户输入
        new_state = dict(state)
        new_state["messages"] = state.get("messages", []) + [HumanMessage(content=user_input)]
        
        return {"messages": [HumanMessage(content=user_input)]}
    
    except KeyboardInterrupt:
        logger.info("User interrupted the session")
        return {"messages": [HumanMessage(content="exit")]}
    except Exception as e:
        logger.error(f"Error getting user input: {str(e)}")
        return {"messages": [HumanMessage(content="发生错误，请重试")]}

def normal_mode_node(state: AgentState) -> Dict[str, Any]:
    """
    Normal模式主节点函数
    
    应用组合模式：组合多个处理器完成复杂任务
    """
    logger.info("Entering normal mode node")
    
    try:
        # 创建处理器实例
        handler = NormalModeHandler()
        
        # 如果没有消息，提供初始指导
        if not state.get("messages"):
            logger.info("No messages found, providing initial guidance")
            return handler.provide_guidance(state)
        
        # 检查是否需要处理模式切换
        mode_switch_result = handler.handle_mode_switch_request(state)
        if mode_switch_result:
            logger.info("Mode switch detected")
            return mode_switch_result
        
        # 处理正常的用户输入
        result = handler.process_user_input(state)
        
        logger.info("Normal mode processing completed")
        return result
    
    except Exception as e:
        logger.error(f"Error in normal mode node: {str(e)}")
        error_message = AIMessage(
            content="抱歉，系统出现错误。请重试或重新启动程序。"
        )
        return {"messages": [error_message]}

def should_continue_in_normal_mode(state: AgentState) -> bool:
    """
    判断是否应该继续在normal模式
    
    应用KISS原则：简单的模式判断逻辑
    """
    current_mode = state.get("mode", "normal")
    
    # 检查是否有模式切换的工具调用
    if state.get("messages"):
        last_message = state["messages"][-1]
        if hasattr(last_message, "tool_calls") and last_message.tool_calls:
            for tool_call in last_message.tool_calls:
                if tool_call.get("name") in ["switch_to_plan_mode", "switch_to_execute_mode"]:
                    return False
    
    return current_mode == "normal"

def handle_normal_mode_tools(state: AgentState) -> Dict[str, Any]:
    """
    处理normal模式下的工具调用结果
    
    遵循DRY原则：统一的工具调用结果处理
    """
    try:
        # 这里可以添加特定于normal模式的工具调用后处理逻辑
        # 目前保持简单，直接返回状态
        logger.info("Processing normal mode tool results")
        
        # 检查工具调用结果，如果有错误需要特殊处理
        if state.get("messages"):
            last_message = state["messages"][-1]
            if hasattr(last_message, "content") and "错误" in str(last_message.content):
                logger.warning("Tool call resulted in error")
        
        return {}
    
    except Exception as e:
        logger.error(f"Error handling normal mode tools: {str(e)}")
        return {}

# ============================================================================
# 便捷函数和工具 - 遵循DRY原则
# ============================================================================

def create_welcome_message(validation_results=None) -> AIMessage:
    """
    创建欢迎消息
    
    应用工厂模式：统一的消息创建
    """
    welcome_content = """🧬 **RNA-seq智能助手** 已启动！

我是您的专业RNA-seq智能助手，可以帮助您：
- 📁 查看和管理FASTQ文件
- 🧬 配置基因组参考文件  
- 📋 制定个性化处理计划
- 🚀 执行完整的RNA-seq流程

请告诉我您的需求，或输入"帮助"查看详细功能介绍。"""

    # 如果有验证结果，添加系统状态信息
    if validation_results:
        status_content = "\n\n📋 **系统验证状态**:\n"
        for status, message in validation_results:
            status_content += f"  {status} {message}\n"
        
        # 统计结果
        success_count = sum(1 for status, _ in validation_results if status == "✅")
        total_count = len(validation_results)
        status_content += f"\n📊 总结: {success_count}/{total_count} 项验证通过"
        
        if success_count == total_count:
            status_content += " 🎉"
        
        welcome_content += status_content
    
    return AIMessage(content=welcome_content)

def create_help_message() -> AIMessage:
    """
    创建帮助消息
    
    应用工厂模式：统一的帮助信息创建
    """
    return AIMessage(content="""
📖 **功能说明**：

🔍 **信息查询命令**：
- "查看目录 [路径]" - 查看指定目录内容
- "查询FASTQ文件 [路径]" - 查看FASTQ文件信息
- "查询基因组 [名称]" - 查看基因组配置（如hg38、mm39）
- "当前配置" - 查看当前nextflow配置

💬 **对话交互**：
- 直接描述您的处理需求
- 询问RNA-seq相关问题
- 寻求处理建议和指导

🚀 **开始处理**：
- "开始分析" - 进入计划制定模式
- "制定计划" - 开始制定处理计划

⚙️ **系统命令**：
- "帮助" - 显示此帮助信息
- "退出" - 退出程序

有什么问题尽管问我！
    """)

def is_help_request(message_content: str) -> bool:
    """
    判断是否为帮助请求
    
    应用KISS原则：简单的帮助请求识别
    """
    help_keywords = ["帮助", "help", "?", "？", "指导", "说明"]
    content_lower = message_content.lower().strip()
    return any(keyword in content_lower for keyword in help_keywords)

def is_analysis_start_request(message_content: str) -> bool:
    """
    判断是否为开始分析请求
    
    应用KISS原则：简单的分析开始请求识别
    """
    start_keywords = ["开始分析", "制定计划", "开始", "分析", "start analysis", "begin"]
    content_lower = message_content.lower().strip()
    return any(keyword in content_lower for keyword in start_keywords)

def _clean_unicode_content(content: str) -> str:
    """
    清理Unicode内容中的无效字符
    
    应用KISS原则：简单有效的字符清理
    """
    try:
        import re
        # 移除代理对字符和其他无效Unicode字符
        cleaned = content.encode('utf-8', errors='ignore').decode('utf-8')
        
        # 进一步清理：移除控制字符但保留换行符和制表符
        cleaned = re.sub(r'[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\x9F]', '', cleaned)
        
        return cleaned
    except Exception as e:
        logger.error(f"Error cleaning unicode content: {str(e)}")
        return "内容包含无效字符，已清理。请重新提供您的需求。"