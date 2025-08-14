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
        处理LLM的结构化响应（JsonOutputParser已处理JSON格式）
        
        返回: (AIMessage, tool_calls列表)
        """
        try:
            # JsonOutputParser已经返回dict格式，直接使用
            if isinstance(response, dict):
                logger.info(f"Plan模式收到结构化响应: {list(response.keys())}")
                
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
            else:
                # 降级处理：如果不是dict，可能是旧格式
                logger.warning(f"Plan模式收到非结构化响应: {type(response)}")
                content = str(response) if response else "响应为空"
                return AIMessage(content=content), []
            
        except Exception as e:
            logger.error(f"Plan模式处理响应时出错: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            return AIMessage(content="处理响应时出现错误"), []
    
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
            
            # 构建计划请求消息
            plan_request = f"""
用户请求制定RNA-seq分析计划。用户输入：{user_input}

请严格按照以下步骤制定分析计划：
1. 首先调用 generate_analysis_task_list 工具自动检测本地文件并生成推荐配置
2. **关键步骤**：根据generate_analysis_task_list的结果，必须调用update_nextflow_param或batch_update_nextflow_config工具实际保存配置到系统状态
3. 调用get_current_nextflow_config验证配置已正确保存
4. 基于最终配置制定详细的分析计划
5. 向用户展示完整的执行流程和当前配置状态
6. 询问用户是否确认或需要修改

⚠️ **重要要求**：
- 必须先调用 generate_analysis_task_list 工具获取智能配置建议
- **必须**调用 update_nextflow_param 或 batch_update_nextflow_config 实际保存配置
- 不能只在回复中显示配置，必须实际更新系统状态
- 每个检测到的配置项都要调用工具保存
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
                "plan_status": "draft"
            }
        
        except Exception as e:
            logger.error(f"Error processing plan request: {str(e)}")
            error_message = AIMessage(
                content=f"制定计划时出现错误：{str(e)}。请重试或提供更多信息。"
            )
            return {"messages": [error_message]}
    
    def handle_plan_modification(self, state: AgentState) -> Dict[str, Any]:
        """
        处理计划修改请求，使用JSON-first架构
        
        遵循开放封闭原则：易于扩展新的修改类型
        """
        try:
            # 获取用户输入
            user_input = ""
            if state.get("messages"):
                last_message = state["messages"][-1]
                if hasattr(last_message, "content"):
                    user_input = last_message.content
            
            # 调用LLM处理修改请求
            response = self.chain.invoke({
                "messages": state["messages"],
                "input": user_input
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
    Plan模式主节点函数，采用JSON-first架构
    
    应用策略模式：根据不同情况采用不同处理策略
    """
    logger.info("Entering plan mode node with JSON-first architecture")
    
    try:
        # 获取UI管理器
        ui_manager = get_ui_manager()
        
        # 显示模式切换信息
        if state.get("mode") != "plan":
            ui_manager.show_mode_switch(state.get("mode", "normal"), "plan", "开始制定分析计划")
        
        # 创建处理器实例
        handler = PlanModeHandler()
        
        # 检查是否需要处理模式切换
        mode_switch_result = handler.handle_mode_switch_request(state)
        if mode_switch_result:
            logger.info("Mode switch detected in plan mode")
            if mode_switch_result.get("mode") == "execute":
                ui_manager.show_mode_switch("plan", "execute", "计划已确认，开始执行")
            return mode_switch_result
        
        # 检查是否是初次进入plan模式（需要制定新计划）
        current_plan = state.get("plan", [])
        plan_status = state.get("plan_status", "")
        
        # 只有在没有计划或明确要求重新制定时才创建新计划
        if not current_plan and plan_status != "created":
            logger.info("Creating new analysis plan using intelligent task list")
            ui_manager.show_info("正在制定RNA-seq分析计划...")
            result = handler.process_plan_request(state)
            # 标记计划已创建，避免重复制定
            result["plan_status"] = "created"
            result["mode"] = "plan"  # 确保模式正确
            return result
        else:
            logger.info("Plan already exists, handling user input with JSON architecture")
            # 处理用户输入（修改计划或确认执行）
            result = handler.handle_plan_modification(state)
            result["mode"] = "plan"
            return result
    
    except Exception as e:
        logger.error(f"Error in plan mode node: {str(e)}")
        error_message = AIMessage(
            content="计划模式出现错误。请重试或返回normal模式。"
        )
        return {
            "messages": [error_message],
            "mode": "plan"
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