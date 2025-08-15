"""
Plan Mode节点 - 制定分析计划和修改nextflow参数
遵循单一职责原则：专门处理plan模式下的计划制定和参数配置
采用JSON-first架构，与normal模式保持一致
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
        # 使用结构化链用于JSON格式输出
        self.chain = create_structured_chain_for_mode("plan")
    
    def _process_llm_response(self, response) -> Tuple[AIMessage, List[Dict[str, Any]]]:
        """
        处理LLM的结构化响应（.with_structured_output()返回Pydantic模型实例）
        增强的错误处理和文本响应转换
        
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
                    logger.info(f"Plan模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            
            # 处理情况2：字符串响应（Gemini有时返回文本而不是结构化格式）
            elif isinstance(response, str):
                logger.warning(f"Plan模式收到意外的字符串响应，尝试智能解析")
                
                # 尝试从文本中提取有用信息
                try:
                    # 检查是否包含JSON格式的配置信息
                    import re
                    import json
                    
                    # 搜索JSON配置
                    json_matches = re.findall(r'\{[^}]*\}', response)
                    config_changes = {}
                    
                    for match in json_matches:
                        try:
                            parsed = json.loads(match)
                            if isinstance(parsed, dict):
                                config_changes.update(parsed)
                        except:
                            continue
                    
                    # 创建基本的结构化响应
                    ai_message = AIMessage(content=response)
                    
                    # 如果找到配置更改，创建工具调用
                    tool_calls = []
                    if config_changes:
                        logger.info(f"从文本响应中提取到配置更改: {config_changes}")
                        # 这里可以根据需要添加工具调用逻辑
                    
                    return ai_message, tool_calls
                    
                except Exception as e:
                    logger.warning(f"解析文本响应失败: {e}")
                    # 降级处理：直接返回文本内容
                    return AIMessage(content=response), []
            
            # 处理情况3：兼容旧的dict格式返回
            elif isinstance(response, dict):
                logger.info(f"Plan模式收到dict格式响应: {list(response.keys())}")
                
                # 提取响应信息
                user_message = response.get("response", "计划制定完成")
                plan_steps = response.get("plan_steps", [])
                config_changes = response.get("config_changes", {})
                ready_to_execute = response.get("ready_to_execute", False)
                tool_calls = response.get("tool_calls", [])
                
                # 构建详细响应
                detailed_response = user_message
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
                        tool_call_obj = {
                            "name": tool_call.get("tool_name"),
                            "args": tool_call.get("parameters", {}),
                            "id": f"call_plan_{i}",
                            "type": "tool_call"
                        }
                        langchain_tool_calls.append(tool_call_obj)
                    
                    ai_message.tool_calls = langchain_tool_calls
                    logger.info(f"Plan模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            
            # 处理情况4：其他类型的响应
            else:
                # 降级处理：如果不是期望的格式
                logger.warning(f"Plan模式收到未知响应格式: {type(response)}")
                content = str(response) if response else "响应为空"
                return AIMessage(content=content), []
            
        except Exception as e:
            logger.error(f"Plan模式处理响应时出错: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            
            # 错误恢复：尝试从原始响应中提取文本内容
            try:
                if hasattr(response, 'content'):
                    content = response.content
                elif isinstance(response, str):
                    content = response
                else:
                    content = str(response)
                
                return AIMessage(content=f"处理响应时出现错误，原始内容: {content}"), []
            except:
                return AIMessage(content=f"处理响应时出现严重错误: {str(e)}"), []
    
    def analyze_requirements(self, state: AgentState) -> Dict[str, Any]:
        """
        分析用户需求和现有信息
        
        应用KISS原则：简单直接的需求分析
        """
        try:
            # 收集当前状态信息
            fastq_info = state.get("fastq_info", {})
            genome_info = state.get("genome_info", {})
            current_config = state.get("nextflow_config", {})
            
            analysis_summary = {
                "has_fastq_info": bool(fastq_info),
                "has_genome_info": bool(genome_info),
                "current_config": current_config,
                "missing_info": []
            }
            
            # 检查缺失的关键信息
            if not fastq_info:
                analysis_summary["missing_info"].append("FASTQ文件信息")
            if not genome_info:
                analysis_summary["missing_info"].append("基因组配置信息")
            
            logger.info(f"Requirements analysis: {analysis_summary}")
            return analysis_summary
        
        except Exception as e:
            logger.error(f"Error analyzing requirements: {str(e)}")
            return {"error": str(e)}
    
    def create_analysis_plan(self, state: AgentState, requirements: Dict[str, Any]) -> List[str]:
        """
        创建分析计划
        
        遵循DRY原则：基于模板的计划生成
        """
        try:
            plan_steps = []
            
            # 基础步骤模板
            base_steps = [
                "1. 数据质量控制 (FastP)",
                "2. 序列比对 (STAR)",
                "3. 基因定量 (featureCounts)",
                "4. 结果整理和报告生成"
            ]
            
            # 根据需求添加前置步骤
            if "FASTQ文件信息" in requirements.get("missing_info", []):
                plan_steps.append("0. 收集和验证FASTQ文件信息")
            
            if "基因组配置信息" in requirements.get("missing_info", []):
                plan_steps.append("0. 配置基因组参考文件")
            
            # 添加基础分析步骤
            plan_steps.extend(base_steps)
            
            logger.info(f"Created analysis plan with {len(plan_steps)} steps")
            return plan_steps
        
        except Exception as e:
            logger.error(f"Error creating analysis plan: {str(e)}")
            return ["错误：无法创建分析计划"]
    
    def configure_nextflow_params(self, state: AgentState, plan_steps: List[str]) -> Dict[str, Any]:
        """
        配置nextflow参数
        
        应用SOLID原则：基于计划步骤配置相应参数
        """
        try:
            config_updates = {}
            
            # 基于计划步骤确定需要启用的流程
            if any("FastP" in step for step in plan_steps):
                config_updates["run_fastp"] = True
            
            if any("STAR" in step for step in plan_steps):
                config_updates["run_star_align"] = True
                config_updates["run_build_star_index"] = True
            
            if any("featureCounts" in step for step in plan_steps):
                config_updates["run_featurecounts"] = True
            
            # 设置输出目录
            config_updates["data"] = "./data"
            
            logger.info(f"Configured nextflow parameters: {config_updates}")
            return config_updates
        
        except Exception as e:
            logger.error(f"Error configuring nextflow params: {str(e)}")
            return {}
    
    def process_plan_request(self, state: AgentState) -> Dict[str, Any]:
        """
        处理计划制定请求，使用智能任务列表生成
        
        应用JSON-first架构：结构化的计划制定流程
        """
        try:
            # 确保当前处于plan模式
            if state.get("mode") != "plan":
                logger.warning(f"Expected plan mode, but got {state.get('mode')}")
                state = update_state_mode(state, "plan")
            
            # 获取最后一条用户消息
            user_input = ""
            if state.get("messages"):
                last_message = state["messages"][-1]
                if hasattr(last_message, "content"):
                    user_input = last_message.content
            
            # 构建初始计划请求消息
            plan_request = f"""
用户请求制定RNA-seq分析计划。用户输入：{user_input}

请调用 generate_analysis_task_list 工具自动检测本地文件并生成推荐配置。
这是制定计划的第一步，稍后我会基于检测结果为用户生成完整的分析计划。
            """
            
            # 调用LLM处理计划请求
            response = self.chain.invoke({
                "messages": state["messages"] + [HumanMessage(content=plan_request)]
            })
            
            # 处理LLM的结构化响应
            parsed_response, _ = self._process_llm_response(response)
            
            logger.info("Plan request processed successfully with JSON response")
            logger.info(f"Plan模式返回的消息tool_calls属性: {hasattr(parsed_response, 'tool_calls')} - {getattr(parsed_response, 'tool_calls', None)}")
            
            # 返回结果，包含解析后的响应
            return {
                "messages": [parsed_response],
                "mode": "plan",
                "plan_status": "initial_request"  # 标记为初始请求，待工具执行
            }
        
        except Exception as e:
            logger.error(f"Error processing plan request: {str(e)}")
            error_message = AIMessage(
                content=f"制定计划时出现错误：{str(e)}。请重试或提供更多信息。"
            )
            return {"messages": [error_message]}
    
    def process_tool_results(self, state: AgentState) -> Dict[str, Any]:
        """
        处理工具执行结果，生成用户友好的计划总结
        
        这是工具调用后的第二轮LLM处理环节
        """
        try:
            logger.info("处理工具执行结果，生成计划总结")
            
            # 检查工具执行结果
            if not state.get("messages"):
                logger.warning("没有工具执行结果")
                return {"messages": [AIMessage(content="没有检测到工具执行结果")]}
            
            # 获取当前配置状态
            current_config = state.get("nextflow_config", {})
            fastq_info = state.get("fastq_info", {})
            genome_info = state.get("genome_info", {})
            
            # 构建工具结果总结请求
            tool_result_request = f"""
工具执行已完成，现在请基于以下信息为用户生成完整的RNA-seq分析计划：

**当前配置状态：**
{current_config}

**FASTQ文件信息：**
{fastq_info}

**基因组信息：**
{genome_info}

请生成一个用户友好的分析计划总结，包括：
1. 检测到的文件和配置概况
2. 推荐的分析流程步骤
3. 当前的关键配置参数
4. 是否需要用户进一步确认或修改

请直接提供格式化的用户回复，不需要再调用工具。
            """
            
            # 调用LLM生成最终的计划总结
            response = self.chain.invoke({
                "messages": state["messages"] + [HumanMessage(content=tool_result_request)]
            })
            
            # 处理LLM的响应
            parsed_response, _ = self._process_llm_response(response)
            
            logger.info("Tool results processed successfully")
            
            return {
                "messages": [parsed_response],
                "mode": "plan",
                "plan_status": "completed"  # 标记计划已完成
            }
            
        except Exception as e:
            logger.error(f"Error processing tool results: {str(e)}")
            error_message = AIMessage(
                content=f"处理工具结果时出现错误：{str(e)}。请重试。"
            )
            return {"messages": [error_message]}
    
    def handle_plan_modification(self, state: AgentState) -> Dict[str, Any]:
        """
        处理计划修改请求，使用JSON-first架构
        
        遵循开放封闭原则：易于扩展新的修改类型
        """
        try:
            # 调用LLM处理修改请求
            response = self.chain.invoke({
                "messages": state["messages"]
            })
            
            # 处理LLM的结构化响应
            parsed_response, _ = self._process_llm_response(response)
            
            logger.info("Plan modification handled with JSON response")
            return {"messages": [parsed_response]}
        
        except Exception as e:
            logger.error(f"Error handling plan modification: {str(e)}")
            error_message = AIMessage(
                content=f"修改计划时出现错误：{str(e)}。请重试。"
            )
            return {"messages": [error_message]}
    
    def generate_task_list(self, state: AgentState, planning_details: Dict[str, Any]) -> Any:
        """
        调用LLM生成任务列表
        
        基于当前状态和规划详情，让LLM生成下一轮需要执行的任务
        """
        try:
            # 构建包含丰富上下文的提示
            context_info = self._build_planning_context(state, planning_details)
            
            planning_request = f"""
基于当前状态分析，请制定下一轮需要执行的任务。

**当前规划状态：**
- 规划轮数: {planning_details.get('iteration_count', 0)}/5
- 配置完整度: {planning_details.get('config_completeness', 0)}%
- 已执行工具: {planning_details.get('tools_executed', [])}
- 最近决策: {planning_details.get('last_decisions', [])}

**当前配置状态：**
{context_info}

**重要提醒：**
- 不要重复执行已完成的工具: {planning_details.get('tools_executed', [])}
- 如果配置完整度已足够或已达到最大轮数，不要生成新任务
- 优先检测FASTQ文件和基因组信息
- 每轮只生成1-3个最关键的任务

请调用必要的工具来收集缺失信息或完善配置。
            """
            
            # 调用LLM处理规划请求  
            response = self.chain.invoke({
                "messages": state["messages"] + [HumanMessage(content=planning_request)]
            })
            
            logger.info("LLM生成任务列表成功")
            return response
        
        except Exception as e:
            logger.error(f"生成任务列表时出错: {str(e)}")
            return None
    
    def create_final_plan_summary(self, state: AgentState, completion_summary: Dict[str, Any]) -> Dict[str, Any]:
        """
        创建最终的计划总结，展示给用户
        """
        try:
            # 构建详细的总结信息
            summary_request = f"""
规划周期已完成，请为用户生成一个详细的RNA-seq分析计划总结。

**规划完成信息：**
- 总轮数: {completion_summary.get('total_iterations', 0)}
- 最终完整度: {completion_summary.get('final_completeness', 0)}%
- 执行任务数: {completion_summary.get('tasks_executed', 0)}
- 使用工具: {completion_summary.get('tools_used', [])}
- 完成原因: {completion_summary.get('completion_reason', '')}

**最终配置状态：**
{self._format_config_summary(completion_summary.get('final_config', {}))}

**收集到的信息：**
{self._format_collected_info(completion_summary.get('collected_info', {}))}

请生成一个用户友好的计划总结，包括：
1. 检测到的文件和配置概况  
2. 推荐的分析流程步骤
3. 当前的关键配置参数
4. 是否准备好执行，或还需要用户确认什么

请提供格式化的用户回复，让用户清楚了解当前状态。
            """
            
            # 调用LLM生成最终总结
            response = self.chain.invoke({
                "messages": state["messages"] + [HumanMessage(content=summary_request)]
            })
            
            # 处理LLM的响应
            parsed_response, _ = self._process_llm_response(response)
            
            logger.info("创建最终计划总结成功")
            
            return {
                "messages": [parsed_response],
                "completed_tasks": state.get("completed_tasks", []),
                "completion_summary": completion_summary
            }
        
        except Exception as e:
            logger.error(f"创建最终计划总结时出错: {str(e)}")
            error_message = AIMessage(
                content=f"生成计划总结时出现错误：{str(e)}。请重试。"
            )
            return {"messages": [error_message]}
    
    def _build_planning_context(self, state: AgentState, planning_details: Dict[str, Any]) -> str:
        """
        构建规划上下文信息
        """
        try:
            nextflow_config = state.get("nextflow_config", {})
            fastq_info = state.get("fastq_info", {})
            genome_info = state.get("genome_info", {})
            
            context_parts = []
            
            # FASTQ文件状态
            if fastq_info:
                context_parts.append(f"✅ FASTQ文件: 已检测到 {len(fastq_info.get('files', []))} 个文件")
            else:
                context_parts.append("❌ FASTQ文件: 未检测")
            
            # 基因组状态  
            if genome_info:
                context_parts.append(f"✅ 基因组信息: {genome_info.get('selected_genome', '已配置')}")
            else:
                context_parts.append("❌ 基因组信息: 未配置")
            
            # Nextflow配置关键参数
            key_configs = [
                ("数据源", nextflow_config.get("local_fastq_files") or nextflow_config.get("srr_ids") or "未设置"),
                ("基因组", nextflow_config.get("genome_version") or nextflow_config.get("local_genome_path") or "未设置"),
                ("质量控制", "启用" if nextflow_config.get("run_fastp") else "禁用"),
                ("序列比对", "启用" if nextflow_config.get("run_star_align") else "禁用")
            ]
            
            context_parts.append("\n**关键配置:**")
            for key, value in key_configs:
                context_parts.append(f"- {key}: {value}")
            
            return "\n".join(context_parts)
        
        except Exception as e:
            logger.error(f"构建规划上下文时出错: {str(e)}")
            return "上下文信息获取失败"
    
    def _format_config_summary(self, config: Dict[str, Any]) -> str:
        """
        格式化配置总结
        """
        try:
            summary_parts = []
            
            # 数据源
            data_source = config.get("local_fastq_files") or config.get("srr_ids") or "未设置"
            summary_parts.append(f"- 数据源: {data_source}")
            
            # 基因组
            genome = config.get("genome_version") or config.get("local_genome_path") or "未设置"
            summary_parts.append(f"- 基因组: {genome}")
            
            # 启用的流程
            enabled_processes = []
            process_map = {
                "run_fastp": "质量控制",
                "run_star_align": "序列比对", 
                "run_featurecounts": "基因定量"
            }
            
            for param, name in process_map.items():
                if config.get(param):
                    enabled_processes.append(name)
            
            summary_parts.append(f"- 启用流程: {', '.join(enabled_processes) if enabled_processes else '无'}")
            
            return "\n".join(summary_parts)
        
        except Exception as e:
            logger.error(f"格式化配置总结时出错: {str(e)}")
            return "配置信息获取失败"
    
    def _format_collected_info(self, collected_info: Dict[str, Any]) -> str:
        """
        格式化收集的信息
        """
        try:
            info_parts = []
            
            fastq_info = collected_info.get("fastq_info", {})
            if fastq_info:
                file_count = len(fastq_info.get("files", []))
                info_parts.append(f"- FASTQ文件: 检测到 {file_count} 个文件")
            
            genome_info = collected_info.get("genome_info", {})
            if genome_info:
                available_genomes = len(genome_info.get("available_genomes", []))
                info_parts.append(f"- 基因组: {available_genomes} 个可用基因组")
            
            if not info_parts:
                info_parts.append("- 暂无详细信息")
            
            return "\n".join(info_parts)
        
        except Exception as e:
            logger.error(f"格式化收集信息时出错: {str(e)}")
            return "信息格式化失败"
    
    def handle_mode_switch_request(self, state: AgentState) -> Dict[str, Any]:
        """
        处理模式切换请求，统一使用工具调用检测
        
        遵循单一职责原则：专门处理模式切换，与normal模式保持一致
        """
        try:
            # 检查最后一条消息是否包含模式切换工具调用
            if not state.get("messages"):
                return {}
                
            last_message = state["messages"][-1]
            
            # 只检查工具调用（LLM调用的工具），与normal模式保持一致
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                for tool_call in last_message.tool_calls:
                    if tool_call.get("name") == "switch_to_execute_mode":
                        logger.info("Switching to execute mode requested via tool call")
                        # 保持完整状态，只更新模式
                        result = dict(state)  # 复制现有状态
                        result["mode"] = "execute"
                        result["execution_status"] = "idle"  # 准备执行
                        result["messages"] = state["messages"] + [
                            AIMessage(content="🔄 计划已确认！正在切换到执行模式...")
                        ]
                        return result
            
            return {}
        
        except Exception as e:
            logger.error(f"Error handling mode switch: {str(e)}")
            return {}

def plan_mode_node(state: AgentState) -> Dict[str, Any]:
    """
    Plan模式主节点函数，基于Plan-Execute任务队列系统
    
    实现两阶段流程：
    1. LLM自主的Plan-Execute循环（无用户参与）
    2. 结果展示和用户交互
    """
    from ..state import (
        update_execution_phase, increment_iteration_count, 
        add_decision_to_history, clear_task_queue, update_config_completeness
    )
    from ..task_executor import (
        create_auto_executor, create_completion_judge, create_task_queue_manager
    )
    
    logger.info("Entering plan mode node with Plan-Execute task queue system")
    
    try:
        # 获取UI管理器
        ui_manager = get_ui_manager()
        
        # 获取当前执行阶段
        execution_phase = state.get("execution_phase", "interacting")
        current_mode = state.get("mode", "normal")
        
        logger.info(f"Plan模式当前阶段: {execution_phase}, 当前模式: {current_mode}")
        
        # 创建核心组件
        handler = PlanModeHandler()
        auto_executor = create_auto_executor()
        completion_judge = create_completion_judge()
        task_manager = create_task_queue_manager()
        
        # 检查是否需要处理模式切换（优先级最高）
        mode_switch_result = handler.handle_mode_switch_request(state)
        if mode_switch_result:
            logger.info("Mode switch detected in plan mode")
            if mode_switch_result.get("mode") == "execute":
                ui_manager.show_mode_switch("plan", "execute", "计划已确认，开始执行")
            return mode_switch_result
        
        # === 阶段1：LLM自主的Plan-Execute循环 ===
        
        if execution_phase == "interacting" and current_mode != "plan":
            # 首次进入Plan模式 - 简化流程，直接进行配置初始化和计划生成
            ui_manager.show_mode_switch(current_mode, "plan", "开始制定分析计划")
            
            logger.info("首次进入Plan模式，开始智能配置初始化")
            
            # 直接调用 generate_analysis_task_list 进行配置初始化
            from ..tools import generate_analysis_task_list
            try:
                task_result = generate_analysis_task_list.invoke({"analysis_type": "standard"})
                logger.info("自动配置检测完成")
                ui_manager.show_info("✅ 已完成文件检测和配置初始化")
            except Exception as e:
                logger.error(f"配置初始化失败: {str(e)}")
                task_result = f"配置检测出现问题: {str(e)}"
            
            # 基于配置结果，直接调用计划生成
            plan_response = handler.process_plan_request(state)
            
            # 更新状态为交互阶段，让用户可以查看和修改计划
            plan_response["mode"] = "plan" 
            plan_response["execution_phase"] = "interacting"
            
            logger.info("计划生成完成，返回用户交互")
            return plan_response
        
        if execution_phase == "planning":
            # 规划阶段：LLM生成任务列表
            logger.info("处理规划阶段 - LLM生成任务列表")
            
            # 检查是否应该继续规划
            should_continue, reason, details = completion_judge.should_continue_planning_cycle(state)
            
            if not should_continue:
                # 规划完成，进入结果收集阶段
                logger.info(f"规划完成: {reason}")
                ui_manager.show_info(f"规划完成: {reason}")
                
                state = update_execution_phase(state, "collecting", reason)
                completion_summary = completion_judge.create_completion_summary(state)
                
                # 生成最终用户展示
                result = handler.create_final_plan_summary(state, completion_summary)
                result["mode"] = "plan"
                result["execution_phase"] = "interacting"  # 转入用户交互阶段
                return result
            
            else:
                # 继续规划：增加迭代计数并调用LLM
                state = increment_iteration_count(state)
                ui_manager.show_info(f"第 {details['iteration_count']+1} 轮规划中...")
                
                # 调用LLM生成任务列表
                llm_response = handler.generate_task_list(state, details)
                
                # 解析LLM响应中的工具调用，转换为任务
                if hasattr(llm_response, "tool_calls") and llm_response.tool_calls:
                    tasks = task_manager.create_tasks_from_tool_calls(llm_response.tool_calls)
                    state = task_manager.add_tasks_to_queue(state, tasks)
                    
                    # 记录决策到历史
                    decision = f"生成{len(tasks)}个任务: {[task['tool_name'] for task in tasks]}"
                    state = add_decision_to_history(state, decision)
                    
                    # 转入执行阶段
                    state = update_execution_phase(state, "executing", f"生成{len(tasks)}个任务")
                    execution_phase = "executing"
                else:
                    # LLM没有生成工具调用，可能认为已完成
                    logger.info("LLM未生成工具调用，可能认为规划已完成")
                    state = update_execution_phase(state, "collecting", "LLM未生成新任务")
                    execution_phase = "collecting"
        
        if execution_phase == "executing":
            # 执行阶段：自动执行任务队列
            logger.info("处理执行阶段 - 自动执行任务队列")
            ui_manager.show_info("正在执行任务...")
            
            # 自动执行所有任务
            state = auto_executor.execute_all_tasks(state)
            
            # 执行完成后检查配置完整度
            state = update_config_completeness(state)
            
            # 转回规划阶段继续下一轮（或结束）
            state = update_execution_phase(state, "planning", "任务执行完成，检查是否需要继续规划")
            execution_phase = "planning"
        
        if execution_phase == "collecting":
            # 结果收集阶段：整理所有执行结果
            logger.info("处理结果收集阶段")
            ui_manager.show_info("正在整理分析结果...")
            
            completion_summary = completion_judge.create_completion_summary(state)
            result = handler.create_final_plan_summary(state, completion_summary)
            result["mode"] = "plan"
            result["execution_phase"] = "interacting"  # 转入用户交互阶段
            return result
        
        # === 阶段2：结果展示和用户交互 ===
        
        if execution_phase == "interacting":
            # 用户交互阶段：处理用户反馈
            logger.info("处理用户交互阶段")
            
            # 检查用户是否要求重新规划
            if state.get("messages"):
                last_message = state["messages"][-1]
                if hasattr(last_message, "content"):
                    content = last_message.content.lower()
                    if any(keyword in content for keyword in ["重新规划", "重新制定", "再次规划"]):
                        # 用户要求重新规划，清空队列并重新开始
                        logger.info("用户要求重新规划")
                        state = clear_task_queue(state)
                        state = update_execution_phase(state, "planning", "用户要求重新规划")
                        execution_phase = "planning"
                        # 继续处理规划阶段
                        return plan_mode_node(state)
            
            # 处理用户的计划修改请求
            result = handler.handle_plan_modification(state)
            result["mode"] = "plan"
            result["execution_phase"] = "interacting"
            return result
        
        # 默认情况（不应该到达这里）
        logger.warning(f"未知执行阶段: {execution_phase}")
        return {"mode": "plan", "execution_phase": "interacting"}
    
    except Exception as e:
        logger.error(f"Error in plan mode node: {str(e)}")
        import traceback
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        
        error_message = AIMessage(
            content="计划模式出现错误。请重试或返回normal模式。"
        )
        return {
            "messages": [error_message],
            "mode": "plan",
            "execution_phase": "interacting"  # 出错时重置到交互阶段
        }

def validate_plan_completeness(state: AgentState) -> bool:
    """
    验证计划完整性
    
    应用KISS原则：简单的完整性检查
    """
    try:
        plan = state.get("plan", [])
        nextflow_config = state.get("nextflow_config", {})
        
        # 检查基本要素
        has_plan_steps = len(plan) > 0
        has_fastq_config = bool(nextflow_config.get("local_fastq_files") or 
                               nextflow_config.get("srr_ids"))
        has_genome_config = bool(nextflow_config.get("local_genome_path") or 
                                nextflow_config.get("download_genome_url"))
        has_enabled_processes = any(nextflow_config.get(key, False) 
                                  for key in ["run_fastp", "run_star_align", "run_featurecounts"])
        
        is_complete = has_plan_steps and has_fastq_config and has_genome_config and has_enabled_processes
        
        logger.info(f"Plan completeness check: {is_complete}")
        return is_complete
    
    except Exception as e:
        logger.error(f"Error validating plan completeness: {str(e)}")
        return False

def should_continue_in_plan_mode(state: AgentState) -> bool:
    """
    判断是否应该继续在plan模式
    
    应用KISS原则：简单的模式判断逻辑
    """
    current_mode = state.get("mode", "normal")
    plan_status = state.get("plan_status", "draft")
    
    # 检查是否有模式切换的工具调用
    if state.get("messages"):
        last_message = state["messages"][-1]
        if hasattr(last_message, "tool_calls") and last_message.tool_calls:
            for tool_call in last_message.tool_calls:
                if tool_call.get("name") == "switch_to_execute_mode":
                    return False
    
    return current_mode == "plan" and plan_status != "confirmed"

def create_plan_summary(state: AgentState) -> AIMessage:
    """
    创建计划总结
    
    应用工厂模式：统一的总结消息创建
    """
    try:
        plan = state.get("plan", [])
        nextflow_config = state.get("nextflow_config", {})
        
        summary_content = "📋 **RNA-seq分析计划总结**\n\n"
        
        # 计划步骤
        if plan:
            summary_content += "**分析步骤：**\n"
            for step in plan:
                summary_content += f"- {step}\n"
            summary_content += "\n"
        
        # 配置信息
        summary_content += "**配置参数：**\n"
        key_configs = {
            "数据目录": nextflow_config.get("data", "未设置"),
            "质量控制": "启用" if nextflow_config.get("run_fastp") else "禁用",
            "序列比对": "启用" if nextflow_config.get("run_star_align") else "禁用",
            "基因定量": "启用" if nextflow_config.get("run_featurecounts") else "禁用"
        }
        
        for key, value in key_configs.items():
            summary_content += f"- {key}: {value}\n"
        
        summary_content += "\n如需修改计划，请告诉我具体要调整的内容。\n"
        summary_content += "确认无误后，请说\"开始执行\"进入执行模式。"
        
        return AIMessage(content=summary_content)
    
    except Exception as e:
        logger.error(f"Error creating plan summary: {str(e)}")
        return AIMessage(content="无法生成计划总结，请重试。")

# ============================================================================
# 计划模板和配置 - 遵循配置分离原则
# ============================================================================

class PlanTemplate:
    """
    计划模板类
    
    遵循模板方法模式：提供标准的计划模板
    """
    
    @staticmethod
    def get_standard_rnaseq_plan() -> List[str]:
        """标准RNA-seq分析计划"""
        return [
            "1. 数据预处理和质量控制 (FastP)",
            "2. 构建基因组索引 (STAR index)",
            "3. 序列比对到参考基因组 (STAR align)",
            "4. 基因表达定量 (featureCounts)",
            "5. 结果整理和质量报告生成"
        ]
    
    @staticmethod
    def get_minimal_rnaseq_plan() -> List[str]:
        """最小RNA-seq分析计划"""
        return [
            "1. 序列比对 (STAR)",
            "2. 基因定量 (featureCounts)",
            "3. 结果输出"
        ]
    
    @staticmethod
    def get_comprehensive_rnaseq_plan() -> List[str]:
        """全面RNA-seq分析计划"""
        return [
            "1. 原始数据质量评估",
            "2. 数据预处理和质量控制 (FastP)",
            "3. 清洁数据质量再评估",
            "4. 构建基因组索引 (STAR index)",
            "5. 序列比对到参考基因组 (STAR align)",
            "6. 比对质量评估",
            "7. 基因表达定量 (featureCounts)",
            "8. 定量结果质量控制",
            "9. 生成综合分析报告"
        ]