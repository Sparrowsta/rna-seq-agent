"""
Plan Mode节点 - 制定分析计划和修改nextflow参数
遵循单一职责原则：专门处理plan模式下的计划制定和参数配置
"""

import logging
from typing import Dict, Any, List
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_state_mode, update_nextflow_config, add_plan_step
from ..core import create_chain_for_mode, create_structured_chain_for_mode
from ..ui_manager import get_ui_manager

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PlanModeHandler:
    """
    Plan模式处理器
    
    遵循单一职责原则：专门处理plan模式的业务逻辑
    """
    
    def __init__(self):
        self.chain = create_chain_for_mode("plan")
        self.structured_chain = create_structured_chain_for_mode("plan")
    
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
        处理计划制定请求
        
        应用组合模式：组合多个处理步骤
        """
        try:
            # 确保当前处于plan模式
            if state.get("mode") != "plan":
                logger.warning(f"Expected plan mode, but got {state.get('mode')}")
                state = update_state_mode(state, "plan")
            
            # 分析需求
            requirements = self.analyze_requirements(state)
            if "error" in requirements:
                return {"messages": [AIMessage(content=f"分析需求时出错：{requirements['error']}")]}
            
            # 创建分析计划
            plan_steps = self.create_analysis_plan(state, requirements)
            
            # 配置nextflow参数
            config_updates = self.configure_nextflow_params(state, plan_steps)
            
            # 调用LLM生成详细的计划说明
            plan_context = {
                "requirements": requirements,
                "plan_steps": plan_steps,
                "config_updates": config_updates
            }
            
            # 创建包含计划信息的消息
            plan_message = HumanMessage(content=f"""
请基于以下信息制定详细的RNA-seq分析计划：

需求分析：{requirements}
计划步骤：{plan_steps}
配置更新：{config_updates}

请提供详细的计划说明和下一步建议。
            """)
            
            # 临时添加计划消息到状态中
            temp_state = state.copy()
            temp_state["messages"] = state["messages"] + [plan_message]
            
            response = self.chain.invoke({"messages": temp_state["messages"]})
            
            # 清理响应内容中的无效字符
            if hasattr(response, 'content') and response.content:
                cleaned_content = self._clean_unicode_content(response.content)
                response.content = cleaned_content
            
            # 更新状态
            result = {
                "messages": [response],
                "plan": plan_steps,
                "plan_status": "draft",
                "nextflow_config": {**state.get("nextflow_config", {}), **config_updates}
            }
            
            logger.info("Plan request processed successfully")
            return result
        
        except Exception as e:
            logger.error(f"Error processing plan request: {str(e)}")
            error_message = AIMessage(
                content=f"制定计划时出现错误：{str(e)}。请重试或提供更多信息。"
            )
            return {"messages": [error_message]}
    
    def handle_plan_modification(self, state: AgentState) -> Dict[str, Any]:
        """
        处理计划修改请求
        
        遵循开放封闭原则：易于扩展新的修改类型
        """
        try:
            # 调用LLM处理修改请求
            response = self.chain.invoke({"messages": state["messages"]})
            
            # 清理响应内容中的无效字符
            if hasattr(response, 'content') and response.content:
                cleaned_content = self._clean_unicode_content(response.content)
                response.content = cleaned_content
            
            logger.info("Plan modification handled")
            return {"messages": [response]}
        
        except Exception as e:
            logger.error(f"Error handling plan modification: {str(e)}")
            error_message = AIMessage(
                content=f"修改计划时出现错误：{str(e)}。请重试。"
            )
            return {"messages": [error_message]}
    
    def _clean_unicode_content(self, content: str) -> str:
        """
        清理Unicode内容中的无效字符
        
        应用KISS原则：简单有效的字符清理
        """
        try:
            # 移除代理对字符和其他无效Unicode字符
            cleaned = content.encode('utf-8', errors='ignore').decode('utf-8')
            
            # 进一步清理：移除控制字符但保留换行符和制表符
            import re
            cleaned = re.sub(r'[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\x9F]', '', cleaned)
            
            return cleaned
        except Exception as e:
            logger.warning(f"Error cleaning unicode content: {str(e)}")
            return "内容包含无效字符，已清理。请重新提供您的需求。"
    
    def handle_mode_switch_request(self, state: AgentState) -> Dict[str, Any]:
        """
        处理模式切换请求
        
        遵循单一职责原则：专门处理模式切换
        """
        try:
            # 检查最后一条消息是否包含模式切换工具调用
            last_message = state["messages"][-1]
            
            if hasattr(last_message, "tool_calls") and last_message.tool_calls:
                for tool_call in last_message.tool_calls:
                    if tool_call.get("name") == "switch_to_execute_mode":
                        logger.info("Switching to execute mode requested")
                        # 更新状态模式和计划状态
                        return {
                            "mode": "execute",
                            "plan_status": "confirmed",
                            "messages": [AIMessage(content="计划已确认，切换到执行模式...")]
                        }
            
            return {}
        
        except Exception as e:
            logger.error(f"Error handling mode switch: {str(e)}")
            return {}

def plan_mode_node(state: AgentState) -> Dict[str, Any]:
    """
    Plan模式主节点函数
    
    应用策略模式：根据不同情况采用不同处理策略
    """
    logger.info("Entering plan mode node")
    
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
            logger.info("Creating new analysis plan")
            ui_manager.show_info("正在制定RNA-seq分析计划...")
            result = handler.process_plan_request(state)
            # 标记计划已创建，避免重复制定
            result["plan_status"] = "created"
            result["mode"] = "plan"  # 确保模式正确
            return result
        else:
            logger.info("Plan already exists, providing plan summary and waiting for user input")
            # 不再调用LLM，直接提供计划总结
            summary = create_plan_summary(state)
            return {
                "messages": [summary],
                "mode": "plan",
                "plan_status": "ready"  # 标记为准备状态，等待用户确认或修改
            }
    
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

def get_plan_template(complexity: str = "standard") -> List[str]:
    """
    获取计划模板
    
    应用工厂模式：根据复杂度返回相应模板
    """
    templates = {
        "minimal": PlanTemplate.get_minimal_rnaseq_plan,
        "standard": PlanTemplate.get_standard_rnaseq_plan,
        "comprehensive": PlanTemplate.get_comprehensive_rnaseq_plan
    }
    
    template_func = templates.get(complexity, PlanTemplate.get_standard_rnaseq_plan)
    return template_func()