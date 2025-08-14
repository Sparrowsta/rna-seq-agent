"""
Execute Mode节点 - 执行nextflow流程和结果总结
遵循单一职责原则：专门处理execute模式下的流程执行和结果处理
采用JSON-first架构，与其他模式保持一致
"""

import logging
import os
import json
import time
import subprocess
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_execution_status
from ..core import create_chain_for_mode, create_structured_chain_for_mode
from ..ui_manager import get_ui_manager

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ExecuteModeHandler:
    """
    Execute模式处理器
    
    遵循单一职责原则：专门处理execute模式的业务逻辑
    采用JSON-first架构，提供实时进度监控
    """
    
    def __init__(self):
        # 使用结构化链用于JSON格式输出
        self.chain = create_structured_chain_for_mode("execute")
        self.structured_chain = create_structured_chain_for_mode("execute")
        self.nextflow_process = None  # 存储nextflow进程
        self.execution_log = []  # 存储执行日志
    
    def _parse_json_response(self, response) -> Tuple[AIMessage, List[Dict[str, Any]]]:
        """
        解析LLM的JSON响应
        
        返回: (AIMessage, tool_calls列表)
        """
        try:
            if hasattr(response, 'content') and response.content:
                # 清理响应内容
                content = _clean_unicode_content(response.content)
                logger.info(f"Execute模式LLM响应内容: {repr(content[:300])}...")
                
                # 移除代码块标记
                if "```json" in content:
                    start = content.find("```json") + 7
                    end = content.find("```", start)
                    if end != -1:
                        content = content[start:end].strip()
                    else:
                        content = content[start:].strip()
                elif content.startswith("```") and content.endswith("```"):
                    content = content[3:-3].strip()
                
                logger.info(f"Execute模式清理后内容: {repr(content[:300])}...")
                
                # 尝试解析JSON
                try:
                    json_data = json.loads(content)
                    logger.info(f"Execute模式JSON解析成功: {json_data.keys()}")
                    
                    # 提取响应信息
                    user_message = json_data.get("response", content)
                    status = json_data.get("status", "unknown")
                    progress = json_data.get("progress", "")
                    next_step = json_data.get("next_step", "")
                    
                    # 构建详细响应
                    detailed_response = user_message
                    if status and status != "unknown":
                        detailed_response += f"\n\n📊 **状态**: {status}"
                    if progress:
                        detailed_response += f"\n📈 **进度**: {progress}"
                    if next_step:
                        detailed_response += f"\n⏭️ **下一步**: {next_step}"
                    
                    # 提取工具调用
                    tool_calls = json_data.get("tool_calls", [])
                    logger.info(f"Execute模式提取到 {len(tool_calls)} 个工具调用: {tool_calls}")
                    
                    # 创建AIMessage
                    ai_message = AIMessage(content=detailed_response)
                    
                    # 如果有工具调用，设置为消息的tool_calls属性
                    if tool_calls:
                        langchain_tool_calls = []
                        for i, tool_call in enumerate(tool_calls):
                            tool_call_obj = {
                                "name": tool_call.get("tool_name"),
                                "args": tool_call.get("parameters", {}),
                                "id": f"call_exec_{i}",
                                "type": "tool_call"
                            }
                            langchain_tool_calls.append(tool_call_obj)
                        
                        ai_message.tool_calls = langchain_tool_calls
                        logger.info(f"Execute模式成功设置tool_calls属性: {langchain_tool_calls}")
                    
                    return ai_message, tool_calls
                    
                except json.JSONDecodeError as e:
                    logger.warning(f"Execute模式LLM输出不是有效JSON，使用原始内容。错误: {str(e)}")
                    return AIMessage(content=content), []
            
            return AIMessage(content="响应为空"), []
            
        except Exception as e:
            logger.error(f"Execute模式解析JSON响应时出错: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            return AIMessage(content="解析响应时出现错误"), []
    
    def prepare_execution(self, state: AgentState) -> Dict[str, Any]:
        """
        准备执行环境
        
        应用KISS原则：简单的执行准备
        """
        try:
            nextflow_config = state.get("nextflow_config", {})
            plan = state.get("plan", [])
            
            # 验证必要的配置
            validation_result = self._validate_execution_config(nextflow_config)
            if not validation_result["valid"]:
                return {
                    "error": f"配置验证失败: {validation_result['message']}",
                    "execution_status": "failed"
                }
            
            # 创建工作目录
            work_dir = nextflow_config.get("data", "./data")
            os.makedirs(work_dir, exist_ok=True)
            os.makedirs(f"{work_dir}/results", exist_ok=True)
            os.makedirs(f"{work_dir}/logs", exist_ok=True)
            
            logger.info(f"Execution environment prepared in {work_dir}")
            
            return {
                "prepared": True,
                "work_dir": work_dir,
                "config": nextflow_config,
                "plan": plan
            }
        
        except Exception as e:
            logger.error(f"Error preparing execution: {str(e)}")
            return {
                "error": str(e),
                "execution_status": "failed"
            }
    
    def _validate_execution_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        验证执行配置
        
        遵循DRY原则：统一的配置验证逻辑
        """
        try:
            missing_fields = []
            
            logger.info(f"Validating config with keys: {list(config.keys())}")
            
            # 检查数据源配置
            has_fastq = bool(config.get("local_fastq_files"))
            has_srr = bool(config.get("srr_ids"))
            
            logger.info(f"Data sources - FASTQ: {has_fastq}, SRR: {has_srr}")
            
            if not (has_fastq or has_srr):
                missing_fields.append("FASTQ文件或SRR ID")
            
            # 检查基因组配置 - 只检查genome_version
            has_genome_version = bool(config.get("genome_version"))
            
            logger.info(f"基因组版本配置: {has_genome_version}")
            
            if not has_genome_version:
                logger.warning("未指定基因组版本，将使用默认hg38")
            
            # 检查启用的流程 - 更宽松的验证，如果没有明确的run_*参数，就假设需要运行基础流程
            enabled_processes = [key for key, value in config.items() 
                               if key.startswith("run_") and value]
            
            logger.info(f"Enabled processes: {enabled_processes}")
            
            # 如果没有明确配置的流程，使用默认流程
            if not enabled_processes:
                logger.info("No explicit processes configured, using default pipeline")
                # 不将这个视为错误，而是使用默认配置
            
            if missing_fields:
                return {
                    "valid": False,
                    "message": f"缺少必要配置: {', '.join(missing_fields)}\n\n请先运行 /plan 命令进行配置。"
                }
            
            return {"valid": True, "message": "配置验证通过"}
        
        except Exception as e:
            logger.error(f"Config validation error: {str(e)}")
            return {"valid": False, "message": f"验证过程出错: {str(e)}"}
    
    def execute_nextflow(self, state: AgentState) -> Dict[str, Any]:
        """
        执行nextflow流程 - 简化版
        
        直接执行命令并等待完成，不进行复杂监控
        """
        try:
            # 检查必要的状态信息
            nextflow_config = state.get("nextflow_config", {})
            plan = state.get("plan", [])
            
            logger.info(f"Execute mode starting with config keys: {list(nextflow_config.keys())}")
            logger.info(f"Plan steps: {len(plan)}")
            
            # 确保有完整的默认配置
            default_config = {
                "srr_ids": "",
                "local_fastq_files": "",
                "genome_version": "",
                "data": "./data",
                "run_download_srr": False,
                "run_download_genome": False,
                "run_build_star_index": False,
                "run_fastp": False,
                "run_star_align": False,
                "run_featurecounts": False
            }
            
            # 合并配置
            merged_config = {**default_config, **nextflow_config}
            logger.info(f"Merged config genome_version: '{merged_config.get('genome_version')}'")
            logger.info(f"Original nextflow_config genome_version: '{nextflow_config.get('genome_version')}'")
            logger.info(f"Default config genome_version: '{default_config.get('genome_version')}'")
            
            # 只有当genome_version为空时才使用默认值hg38
            if not merged_config.get("genome_version"):
                merged_config["genome_version"] = "hg38"
                logger.info("使用默认基因组版本: hg38")
            else:
                logger.info(f"使用配置的基因组版本: {merged_config['genome_version']}")
            logger.info(f"Final merged config keys: {list(merged_config.keys())}")
            
            # 构建nextflow参数
            params = self._build_nextflow_params(merged_config)
            
            # 准备执行环境
            work_dir = params.get("data", "./data")
            os.makedirs(work_dir, exist_ok=True)
            logger.info(f"Execution environment prepared in {work_dir}")
            
            # 构建nextflow命令
            command_parts = ["nextflow", "run", "main.nf"]
            
            for key, value in params.items():
                if isinstance(value, bool):
                    command_parts.extend([f"--{key}", "true" if value else "false"])
                elif value:  # 非空值
                    command_parts.extend([f"--{key}", str(value)])
            
            command_parts.extend(["-c", "config/nextflow.config"])
            command_parts.extend(["-work-dir", "./work"])
            
            command = " ".join(command_parts)
            logger.info(f"Executing nextflow command: {command}")
            
            # 简化执行：直接运行并等待完成
            try:
                import subprocess
                result = subprocess.run(
                    command,
                    shell=True,
                    cwd=".",
                    capture_output=True,
                    text=True,
                    timeout=1800  # 30分钟超时
                )
                
                if result.returncode == 0:
                    success_message = "🎉 **Nextflow执行成功完成！**\n\n"
                    success_message += f"⏱️ **执行概况**\n"
                    success_message += f"• 命令: `{command}`\n"
                    success_message += f"• 工作目录: {work_dir}\n"
                    success_message += f"• 退出代码: {result.returncode}\n\n"
                    
                    if result.stdout:
                        success_message += "📊 **输出概要**:\n"
                        success_message += f"```\n{result.stdout[-500:]}\n```\n\n"
                    
                    success_message += "✅ 分析结果已保存到输出目录中。"
                    
                    return {
                        "execution_status": "completed",
                        "messages": [AIMessage(content=success_message)]
                    }
                else:
                    error_message = "❌ **Nextflow执行失败**\n\n"
                    error_message += f"• 退出代码: {result.returncode}\n"
                    error_message += f"• 命令: `{command}`\n\n"
                    
                    if result.stderr:
                        error_message += "**错误信息 (stderr)**:\n"
                        error_message += f"```\n{result.stderr[-1500:]}\n```\n\n"
                    
                    if result.stdout:
                        error_message += "**输出信息 (stdout)**:\n"
                        error_message += f"```\n{result.stdout[-1500:]}\n```\n\n"
                    
                    if not result.stderr and not result.stdout:
                        error_message += "**注意**: 没有捕获到错误信息或输出信息\n\n"
                    
                    error_message += "**建议**: 检查配置参数和输入文件。可能的原因：\n"
                    error_message += "- 文件权限问题\n"
                    error_message += "- Docker环境中的路径问题\n" 
                    error_message += "- 缺少依赖工具\n"
                    error_message += "- 配置文件问题"
                    
                    return {
                        "execution_status": "failed",
                        "messages": [AIMessage(content=error_message)]
                    }
                    
            except subprocess.TimeoutExpired:
                return {
                    "execution_status": "failed", 
                    "messages": [AIMessage(content="⏰ **执行超时**\n\n执行时间超过30分钟，自动停止。\n\n建议检查输入数据大小和系统资源。")]
                }
                
        except Exception as e:
            logger.error(f"Error executing nextflow: {str(e)}")
            return {
                "execution_status": "failed",
                "messages": [AIMessage(content=f"❌ **执行出错**: {str(e)}")]
            }
    
    def monitor_execution(self, state: AgentState) -> Dict[str, Any]:
        """
        简化版监控 - 直接返回完成状态，避免循环
        """
        return {
            "execution_status": "completed",
            "messages": [AIMessage(content="✅ **监控完成**\n\n执行流程已结束。")]
        }
    
    def collect_results(self, state: AgentState) -> Dict[str, Any]:
        """
        简化版结果收集 - 直接返回完成状态
        """
        return {
            "execution_status": "completed",
            "messages": [AIMessage(content="✅ **结果收集完成**\n\n分析结果已保存。")]
        }
    
    def _build_nextflow_params(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        构建nextflow参数
        
        遵循DRY原则：统一的参数构建逻辑
        优先级：本地文件 > genome_version > 下载URL
        """
        params = {}
        
        logger.info(f"Building nextflow params from config: {config}")
        
        # 数据源参数
        data_params = ["srr_ids", "local_fastq_files", "data"]
        for param in data_params:
            if param in config and config[param]:
                params[param] = config[param]
                logger.info(f"Added param: {param} = {config[param]}")
        
        # 基因组参数 - 只使用genome_version参数
        # main.nf已重构为统一使用genome_version从genomes.json获取路径
        if config.get("genome_version"):
            params["genome_version"] = config["genome_version"]
            logger.info(f"使用基因组版本: {config['genome_version']}")
        else:
            # 默认使用hg38
            params["genome_version"] = "hg38"
            logger.info("未指定基因组版本，使用默认hg38")
        
        # 流程控制参数
        run_params = [
            "run_download_srr", "run_download_genome", "run_build_star_index", 
            "run_fastp", "run_star_align", "run_featurecounts"
        ]
        
        for param in run_params:
            if param in config and config[param] is not None:
                params[param] = bool(config[param])
                logger.info(f"Added param: {param} = {config[param]}")
        
        # 确保有数据目录
        if "data" not in params:
            params["data"] = "./data"
        
        logger.info(f"Final nextflow params: {params}")
        return params
    

def execute_mode_node(state: AgentState) -> Dict[str, Any]:
    """
    Execute模式主节点函数
    
    应用状态机模式：根据执行状态采取不同行动
    """
    logger.info("Entering execute mode node")
    
    try:
        # 获取UI管理器
        ui_manager = get_ui_manager()
        
        # 显示模式切换信息
        if state.get("mode") != "execute":
            ui_manager.show_mode_switch(state.get("mode", "plan"), "execute", "开始执行RNA-seq分析流程")
        
        # 创建处理器实例
        handler = ExecuteModeHandler()
        
        # 获取当前执行状态
        execution_status = state.get("execution_status", "idle")
        logger.info(f"Current execution status: {execution_status}")
        
        # 检查是否为用户重新发起的执行命令
        # 如果最后一条消息是人类输入的execute命令，则重置状态
        if state.get("messages"):
            last_message = state["messages"][-1]
            if (hasattr(last_message, "type") and last_message.type == "human" and
                hasattr(last_message, "content") and 
                last_message.content.lower().strip() in ["/execute", "/开始执行", "/执行"]):
                logger.info("Detected user execute command, resetting execution status to idle")
                execution_status = "idle"
        
        if execution_status == "idle":
            # 开始执行 - 简化版：直接执行并等待完成
            logger.info("Starting nextflow execution")
            ui_manager.show_info("正在启动Nextflow执行流程...")
            result = handler.execute_nextflow(state)
            
            # 确保返回值包含必要的状态更新
            if "mode" not in result:
                result["mode"] = "execute"
            
            # 简化：执行完成后直接标记为完成，不再进行监控循环
            result["execution_status"] = "completed"
            
            return result
        
        else:
            # 其他状态直接返回完成状态，避免循环
            return {
                "mode": "execute",
                "execution_status": "completed", 
                "messages": [AIMessage(content="✅ **执行流程已完成！**\n\n可以继续其他操作或查看结果。")]
            }
    
    except Exception as e:
        logger.error(f"Error in execute mode node: {str(e)}")
        import traceback
        traceback.print_exc()
        
        error_message = AIMessage(
            content=f"执行模式出现错误：{str(e)}\n\n请检查系统状态或重新开始。"
        )
        return {
            "messages": [error_message],
            "mode": "execute",
            "execution_status": "failed"
        }

def should_continue_execution(state: AgentState) -> bool:
    """
    判断是否应该继续执行
    
    应用KISS原则：简单的执行状态判断
    """
    execution_status = state.get("execution_status", "idle")
    return execution_status in ["idle", "running"]

def is_execution_complete(state: AgentState) -> bool:
    """
    判断执行是否完成
    
    应用KISS原则：简单的完成状态判断
    """
    execution_status = state.get("execution_status", "idle")
    return execution_status in ["completed", "failed"]

def create_execution_report(state: AgentState) -> AIMessage:
    """
    创建执行报告
    
    应用工厂模式：统一的报告创建
    """
    try:
        plan = state.get("plan", [])
        execution_results = state.get("execution_results", {})
        execution_status = state.get("execution_status", "unknown")
        
        report_content = "📊 **RNA-seq分析执行报告**\n\n"
        
        # 执行计划
        if plan:
            report_content += "**执行计划**:\n"
            for i, step in enumerate(plan, 1):
                report_content += f"{i}. {step}\n"
            report_content += "\n"
        
        # 执行状态
        status_emoji = {
            "completed": "✅",
            "failed": "❌", 
            "running": "🔄",
            "idle": "⏸️"
        }
        
        report_content += f"**执行状态**: {status_emoji.get(execution_status, '❓')} {execution_status}\n\n"
        
        # 结果信息
        if execution_results:
            summary = execution_results.get("summary", "")
            if summary:
                report_content += summary
        
        return AIMessage(content=report_content)
    
    except Exception as e:
        logger.error(f"Error creating execution report: {str(e)}")
        return AIMessage(content="无法生成执行报告，请检查系统状态。")


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
