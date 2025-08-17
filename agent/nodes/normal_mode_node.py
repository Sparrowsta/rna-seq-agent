"""
Normal Mode节点 - 信息收集和用户交互
遵循单一职责原则：专门处理normal模式下的用户交互和信息收集
"""

import logging
from typing import Dict, Any, List, Tuple
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_state_mode
from ..core import create_chain_for_mode, create_dual_llm_chain_for_mode

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NormalModeHandler:
    """
    Normal模式处理器 - 双LLM架构
    
    遵循LangGraph官方最佳实践：分离工具调用和结构化输出
    """
    
    def __init__(self):
        # 获取双LLM链配置
        self.dual_chain = create_dual_llm_chain_for_mode("normal")
        self.tool_chain = self.dual_chain["tool_chain"]
        self.response_chain = self.dual_chain["response_chain"]
    
    def _execute_tool_chain(self, messages, user_input) -> Tuple[Any, List[Dict[str, Any]]]:
        """
        执行工具调用链 - 双LLM架构第一阶段
        
        返回: (llm_response, tool_execution_results列表)
        """
        try:
            logger.info("执行工具调用链...")
            
            # 调用工具链获取工具调用意图
            tool_response = self.tool_chain.invoke({
                "messages": messages,
                "input": user_input
            })
            
            logger.info(f"工具链响应类型: {type(tool_response)}")
            
            # 检查是否有工具调用
            tool_execution_results = []
            if hasattr(tool_response, 'tool_calls') and tool_response.tool_calls:
                logger.info(f"检测到 {len(tool_response.tool_calls)} 个工具调用")
                
                # 实际执行每个工具调用
                from ..core import ALL_TOOLS
                
                # 创建工具名称到工具对象的映射
                tool_map = {tool.name: tool for tool in ALL_TOOLS}
                
                for i, tool_call in enumerate(tool_response.tool_calls):
                    tool_name = tool_call.get('name', 'unknown')
                    tool_args = tool_call.get('args', {})
                    
                    logger.info(f"执行工具调用 {i+1}: {tool_name} with args: {tool_args}")
                    
                    if tool_name in tool_map:
                        try:
                            # 实际执行工具
                            tool_result = tool_map[tool_name].invoke(tool_args)
                            tool_execution_results.append({
                                'tool_name': tool_name,
                                'args': tool_args,
                                'result': tool_result
                            })
                            logger.info(f"工具 {tool_name} 执行成功")
                        except Exception as tool_error:
                            logger.error(f"工具 {tool_name} 执行失败: {str(tool_error)}")
                            tool_execution_results.append({
                                'tool_name': tool_name,
                                'args': tool_args,
                                'result': f"工具执行失败: {str(tool_error)}"
                            })
                    else:
                        logger.error(f"未知工具: {tool_name}")
                        tool_execution_results.append({
                            'tool_name': tool_name,
                            'args': tool_args,
                            'result': f"未知工具: {tool_name}"
                        })
            else:
                logger.info("没有检测到工具调用")
            
            return tool_response, tool_execution_results
            
        except Exception as e:
            logger.error(f"工具链执行失败: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            raise e
    
    def _execute_response_chain(self, original_input: str, tool_results: str) -> Any:
        """
        执行结构化响应链 - 双LLM架构第二阶段
        
        返回: 结构化响应对象
        """
        try:
            logger.info("执行结构化响应链...")
            
            # 调用结构化响应链
            structured_response = self.response_chain.invoke({
                "original_input": original_input,
                "tool_results": tool_results
            })
            
            logger.info(f"结构化响应类型: {type(structured_response)}")
            logger.info(f"结构化响应内容: {structured_response}")
            
            return structured_response
            
        except Exception as e:
            logger.error(f"结构化响应链执行失败: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            raise e
    
    def _format_tool_results(self, tool_response, tool_execution_results) -> str:
        """
        格式化工具执行结果
        
        将工具执行结果转换为可读的字符串格式，传递给第二阶段LLM
        """
        try:
            if not tool_execution_results:
                return "没有执行任何工具调用。"
            
            formatted_results = []
            
            # 处理实际的工具执行结果
            for result in tool_execution_results:
                tool_name = result.get('tool_name', 'unknown')
                args = result.get('args', {})
                tool_result = result.get('result', '无结果')
                
                formatted_results.append(f"工具: {tool_name}")
                if args:
                    formatted_results.append(f"参数: {args}")
                formatted_results.append(f"执行结果: {tool_result}")
                formatted_results.append("---")
            
            return "\n".join(formatted_results)
            
        except Exception as e:
            logger.error(f"格式化工具结果失败: {str(e)}")
            return f"工具执行完成，但结果格式化出错: {str(e)}"
    
    def _create_final_ai_message(self, structured_response) -> AIMessage:
        """
        基于结构化响应创建最终的AI消息
        """
        try:
            # 检查是否是NormalModeResponse格式
            if hasattr(structured_response, 'response'):
                user_message = getattr(structured_response, "response", "响应完成")
                suggested_actions = getattr(structured_response, "suggested_actions", [])
                need_more_info = getattr(structured_response, "need_more_info", False)
                
                # 构建详细响应
                detailed_response = user_message
                
                if suggested_actions:
                    detailed_response += "\n\n💡 **建议操作：**\n"
                    detailed_response += "\n".join([f"  - {action}" for action in suggested_actions])
                
                if need_more_info:
                    detailed_response += "\n\n❓ 需要更多信息才能继续。"
                
                return AIMessage(content=detailed_response)
            else:
                # 如果不是预期格式，直接使用字符串表示
                return AIMessage(content=str(structured_response))
                
        except Exception as e:
            logger.error(f"创建AI消息失败: {str(e)}")
            return AIMessage(content=f"响应处理出错: {str(e)}")
    
    def process_user_input(self, state: AgentState) -> Dict[str, Any]:
        """
        双阶段处理用户输入 - LangGraph官方双LLM架构
        
        第一阶段：工具调用链处理工具执行
        第二阶段：结构化响应链格式化输出
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
            
            logger.info(f"开始双阶段处理用户输入: {user_input}")
            
            # 第一阶段：执行工具调用链
            try:
                tool_response, tool_calls = self._execute_tool_chain(state["messages"], user_input)
                logger.info(f"第一阶段完成 - 工具调用链执行成功")
            except Exception as tool_error:
                logger.error(f"第一阶段失败 - 工具调用链执行失败: {str(tool_error)}")
                raise tool_error
            
            # 格式化工具执行结果
            tool_results_text = self._format_tool_results(tool_response, tool_calls)
            logger.info(f"工具执行结果格式化完成: {tool_results_text}")
            
            # 第二阶段：执行结构化响应链
            try:
                structured_response = self._execute_response_chain(user_input, tool_results_text)
                logger.info(f"第二阶段完成 - 结构化响应链执行成功")
            except Exception as response_error:
                logger.error(f"第二阶段失败 - 结构化响应链执行失败: {str(response_error)}")
                raise response_error
            
            # 创建最终的AI消息
            final_ai_message = self._create_final_ai_message(structured_response)
            
            logger.info(f"双阶段处理完成，返回最终响应")
            
            return {"messages": [final_ai_message]}
        
        except Exception as e:
            logger.error(f"Error in dual-stage processing: {str(e)}")
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