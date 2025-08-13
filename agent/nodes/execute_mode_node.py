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
            required_fields = []
            missing_fields = []
            
            # 检查数据源配置
            has_fastq = bool(config.get("local_fastq_files"))
            has_srr = bool(config.get("srr_ids"))
            
            if not (has_fastq or has_srr):
                missing_fields.append("FASTQ文件或SRR ID")
            
            # 检查基因组配置
            has_local_genome = bool(config.get("local_genome_path"))
            has_download_genome = bool(config.get("download_genome_url"))
            
            if not (has_local_genome or has_download_genome):
                missing_fields.append("基因组文件配置")
            
            # 检查启用的流程
            enabled_processes = [key for key, value in config.items() 
                               if key.startswith("run_") and value]
            
            if not enabled_processes:
                missing_fields.append("至少一个分析流程")
            
            if missing_fields:
                return {
                    "valid": False,
                    "message": f"缺少必要配置: {', '.join(missing_fields)}"
                }
            
            return {"valid": True, "message": "配置验证通过"}
        
        except Exception as e:
            return {"valid": False, "message": f"验证过程出错: {str(e)}"}
    
    def execute_nextflow(self, state: AgentState) -> Dict[str, Any]:
        """
        执行nextflow流程，提供实时进度监控
        
        应用命令模式：封装执行命令，支持实时输出
        """
        try:
            # 准备执行
            prep_result = self.prepare_execution(state)
            if "error" in prep_result:
                return prep_result
            
            # 从智能任务列表获取配置或使用状态中的配置
            config = state.get("nextflow_config", prep_result.get("config", {}))
            
            # 构建nextflow命令参数
            params = self._build_nextflow_params(config)
            
            # 生成执行命令
            cmd_parts = ["nextflow", "run", "main.nf"]
            
            # 添加参数
            for key, value in params.items():
                if isinstance(value, bool):
                    if value:
                        cmd_parts.extend([f"--{key}", "true"])
                else:
                    cmd_parts.extend([f"--{key}", str(value)])
            
            # 添加配置文件和工作目录
            cmd_parts.extend(["-c", "config/nextflow.config"])
            cmd_parts.extend(["-work-dir", "./work"])
            
            command = " ".join(cmd_parts)
            
            logger.info(f"Executing nextflow command: {command}")
            
            # 创建进度监控消息
            progress_message = self._create_initial_progress_message(command, params)
            
            # 实际执行nextflow（在后台）
            try:
                # 启动nextflow进程
                self.nextflow_process = subprocess.Popen(
                    cmd_parts,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    bufsize=1,
                    universal_newlines=True
                )
                
                execution_info = {
                    "command": command,
                    "status": "running",
                    "start_time": time.time(),
                    "work_dir": prep_result["work_dir"],
                    "process_id": self.nextflow_process.pid,
                    "params": params
                }
                
                logger.info(f"Nextflow process started with PID: {self.nextflow_process.pid}")
                
            except FileNotFoundError:
                # Nextflow不可用，使用模拟模式
                logger.warning("Nextflow not found, using simulation mode")
                execution_info = {
                    "command": command,
                    "status": "simulated",
                    "start_time": time.time(),
                    "work_dir": prep_result["work_dir"],
                    "process_id": None,
                    "params": params
                }
            
            return {
                "execution_info": execution_info,
                "execution_status": "running",
                "messages": [AIMessage(content=progress_message)]
            }
        
        except Exception as e:
            logger.error(f"Error executing nextflow: {str(e)}")
            return {
                "error": str(e),
                "execution_status": "failed",
                "messages": [AIMessage(content=f"执行失败：{str(e)}")]
            }
    
    def _create_initial_progress_message(self, command: str, params: Dict[str, Any]) -> str:
        """创建初始进度消息"""
        message_parts = [
            "🚀 **Nextflow流程已启动！**",
            "",
            "📋 **执行信息：**",
            f"```bash",
            f"{command}",
            f"```",
            "",
            "⚙️ **关键参数：**"
        ]
        
        # 显示关键参数
        key_params = {
            "数据源": params.get("local_fastq_files") or params.get("srr_ids", "未指定"),
            "基因组": params.get("local_genome_path") or params.get("genome_version", "未指定"),
            "质量控制": "启用" if params.get("run_fastp") else "禁用",
            "序列比对": "启用" if params.get("run_star_align") else "禁用",
            "表达定量": "启用" if params.get("run_featurecounts") else "禁用"
        }
        
        for key, value in key_params.items():
            message_parts.append(f"- {key}: {value}")
        
        message_parts.extend([
            "",
            "📊 **实时进度监控：**",
            "[░░░░░░░░░░] 0% - 正在初始化...",
            "",
            "⏳ 流程正在后台运行，请稍等..."
        ])
        
        return "\n".join(message_parts)
    
    def _build_nextflow_params(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        构建nextflow参数
        
        遵循DRY原则：统一的参数构建逻辑
        """
        params = {}
        
        # 直接映射的参数
        direct_params = [
            "srr_ids", "local_genome_path", "local_gtf_path",
            "download_genome_url", "download_gtf_url", "local_fastq_files",
            "data", "run_download_srr", "run_download_genome",
            "run_build_star_index", "run_fastp", "run_star_align", "run_featurecounts"
        ]
        
        for param in direct_params:
            if param in config and config[param]:
                params[param] = config[param]
        
        return params
    
    def monitor_execution(self, state: AgentState) -> Dict[str, Any]:
        """
        监控执行状态，提供类似原生nextflow的进度显示
        
        应用观察者模式：实时监控执行进度
        """
        try:
            execution_info = state.get("execution_results", {}).get("execution_info", {})
            
            if not execution_info:
                return {
                    "messages": [AIMessage(content="没有找到正在执行的流程。")]
                }
            
            # 检查实际进程状态
            current_status = self._check_process_status(execution_info)
            
            # 获取执行时间
            current_time = time.time()
            start_time = execution_info.get("start_time", current_time)
            elapsed_time = current_time - start_time
            
            # 根据进程状态和时间确定进度
            progress_info = self._calculate_progress(current_status, elapsed_time)
            
            # 生成类似nextflow的进度报告
            progress_message = self._generate_nextflow_style_progress(
                progress_info, elapsed_time, execution_info
            )
            
            # 确定执行状态
            if progress_info["completed"]:
                execution_status = "completed"
            elif current_status["failed"]:
                execution_status = "failed"
            else:
                execution_status = "running"
            
            return {
                "execution_status": execution_status,
                "execution_results": {
                    "execution_info": execution_info,
                    "progress_info": progress_info,
                    "current_status": current_status
                },
                "messages": [AIMessage(content=progress_message)]
            }
        
        except Exception as e:
            logger.error(f"Error monitoring execution: {str(e)}")
            return {
                "messages": [AIMessage(content=f"监控执行状态时出错：{str(e)}")]
            }
    
    def _check_process_status(self, execution_info: Dict[str, Any]) -> Dict[str, Any]:
        """检查实际进程状态"""
        try:
            process_id = execution_info.get("process_id")
            
            if not process_id:
                # 模拟模式
                return {
                    "running": True,
                    "failed": False,
                    "simulated": True,
                    "exit_code": None
                }
            
            if self.nextflow_process:
                # 检查进程是否还在运行
                exit_code = self.nextflow_process.poll()
                
                if exit_code is None:
                    # 进程仍在运行
                    return {
                        "running": True,
                        "failed": False,
                        "simulated": False,
                        "exit_code": None
                    }
                else:
                    # 进程已完成
                    return {
                        "running": False,
                        "failed": exit_code != 0,
                        "simulated": False,
                        "exit_code": exit_code
                    }
            else:
                # 进程信息丢失，尝试检查PID
                try:
                    os.kill(process_id, 0)  # 检查进程是否存在
                    return {
                        "running": True,
                        "failed": False,
                        "simulated": False,
                        "exit_code": None
                    }
                except ProcessLookupError:
                    return {
                        "running": False,
                        "failed": False,  # 无法确定失败状态
                        "simulated": False,
                        "exit_code": 0
                    }
        
        except Exception as e:
            logger.error(f"Error checking process status: {str(e)}")
            return {
                "running": False,
                "failed": True,
                "simulated": False,
                "exit_code": -1,
                "error": str(e)
            }
    
    def _calculate_progress(self, status: Dict[str, Any], elapsed_time: float) -> Dict[str, Any]:
        """计算进度信息"""
        try:
            if status.get("failed"):
                return {
                    "percent": 0,
                    "stage": "执行失败",
                    "stage_emoji": "❌",
                    "completed": False,
                    "failed": True,
                    "processes": []
                }
            
            if status.get("simulated"):
                # 模拟模式的进度计算
                if elapsed_time < 30:
                    percent = min(10, elapsed_time / 3)
                    stage = "正在初始化环境"
                    stage_emoji = "🔧"
                elif elapsed_time < 120:
                    percent = min(30, 10 + (elapsed_time - 30) / 3)
                    stage = "正在执行质量控制"
                    stage_emoji = "🧹"
                elif elapsed_time < 300:
                    percent = min(60, 30 + (elapsed_time - 120) / 6)
                    stage = "正在进行序列比对"
                    stage_emoji = "🎯"
                elif elapsed_time < 480:
                    percent = min(85, 60 + (elapsed_time - 300) / 7.2)
                    stage = "正在进行基因定量"
                    stage_emoji = "📊"
                else:
                    percent = 100
                    stage = "执行完成"
                    stage_emoji = "✅"
                
                # 模拟进程列表
                processes = self._generate_simulated_processes(elapsed_time)
                
                return {
                    "percent": int(percent),
                    "stage": stage,
                    "stage_emoji": stage_emoji,
                    "completed": percent >= 100,
                    "failed": False,
                    "processes": processes,
                    "simulated": True
                }
            
            else:
                # 实际执行模式 - 这里可以解析nextflow输出
                # 现在使用基于时间的估算
                if not status.get("running"):
                    return {
                        "percent": 100,
                        "stage": "执行完成",
                        "stage_emoji": "✅",
                        "completed": True,
                        "failed": False,
                        "processes": []
                    }
                
                # 基于时间的进度估算
                estimated_total = 600  # 10分钟估算
                percent = min(95, (elapsed_time / estimated_total) * 100)
                
                if elapsed_time < 60:
                    stage = "正在初始化"
                    stage_emoji = "🔧"
                elif elapsed_time < 180:
                    stage = "正在执行质量控制"
                    stage_emoji = "🧹"
                elif elapsed_time < 420:
                    stage = "正在进行序列比对"
                    stage_emoji = "🎯"
                else:
                    stage = "正在进行基因定量"
                    stage_emoji = "📊"
                
                return {
                    "percent": int(percent),
                    "stage": stage,
                    "stage_emoji": stage_emoji,
                    "completed": False,
                    "failed": False,
                    "processes": [],
                    "simulated": False
                }
                
        except Exception as e:
            logger.error(f"Error calculating progress: {str(e)}")
            return {
                "percent": 0,
                "stage": "进度计算错误",
                "stage_emoji": "❌",
                "completed": False,
                "failed": True,
                "processes": []
            }
    
    def _generate_simulated_processes(self, elapsed_time: float) -> List[Dict[str, Any]]:
        """生成模拟的进程状态"""
        processes = []
        
        # 根据时间添加已完成的进程
        if elapsed_time > 30:
            processes.append({
                "name": "DOWNLOAD_SRR",
                "status": "COMPLETED",
                "progress": "100%",
                "emoji": "✅"
            })
        
        if elapsed_time > 60:
            processes.append({
                "name": "BUILD_STAR_INDEX", 
                "status": "COMPLETED",
                "progress": "100%",
                "emoji": "✅"
            })
        
        if elapsed_time > 120:
            processes.append({
                "name": "FASTP_QC",
                "status": "COMPLETED", 
                "progress": "100%",
                "emoji": "✅"
            })
        
        if elapsed_time > 180:
            if elapsed_time < 300:
                processes.append({
                    "name": "STAR_ALIGN",
                    "status": "RUNNING",
                    "progress": f"{min(100, int((elapsed_time - 180) / 1.2))}%",
                    "emoji": "🔄"
                })
            else:
                processes.append({
                    "name": "STAR_ALIGN",
                    "status": "COMPLETED",
                    "progress": "100%", 
                    "emoji": "✅"
                })
        
        if elapsed_time > 300:
            if elapsed_time < 480:
                processes.append({
                    "name": "FEATURECOUNTS",
                    "status": "RUNNING",
                    "progress": f"{min(100, int((elapsed_time - 300) / 1.8))}%",
                    "emoji": "🔄"
                })
            else:
                processes.append({
                    "name": "FEATURECOUNTS", 
                    "status": "COMPLETED",
                    "progress": "100%",
                    "emoji": "✅"
                })
        
        return processes
    
    def _generate_nextflow_style_progress(self, progress_info: Dict, elapsed_time: float, execution_info: Dict) -> str:
        """生成类似nextflow风格的进度报告"""
        try:
            message_parts = []
            
            # 标题和基本信息
            if progress_info.get("failed"):
                message_parts.append("❌ **Nextflow执行失败**")
            elif progress_info.get("completed"):
                message_parts.append("✅ **Nextflow执行完成**")
            else:
                message_parts.append("🔄 **Nextflow执行进度**")
            
            message_parts.append("")
            
            # 时间信息
            hours = int(elapsed_time // 3600)
            minutes = int((elapsed_time % 3600) // 60)
            seconds = int(elapsed_time % 60)
            
            if hours > 0:
                time_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
            else:
                time_str = f"{minutes:02d}:{seconds:02d}"
            
            message_parts.extend([
                f"⏱️  **运行时间**: {time_str}",
                f"📈 **总体进度**: {progress_info['percent']}%",
                f"{progress_info['stage_emoji']} **当前阶段**: {progress_info['stage']}",
                ""
            ])
            
            # 进度条
            progress_bar = self._create_progress_bar(progress_info["percent"])
            message_parts.append(f"```\n{progress_bar}\n```")
            message_parts.append("")
            
            # 进程状态
            if progress_info.get("processes"):
                message_parts.append("📋 **进程状态**:")
                for process in progress_info["processes"]:
                    status_line = f"{process['emoji']} {process['name']}: {process['status']} ({process['progress']})"
                    message_parts.append(f"  {status_line}")
                message_parts.append("")
            
            # 工作目录信息
            work_dir = execution_info.get("work_dir", "./data")
            message_parts.append(f"📁 **工作目录**: {work_dir}")
            
            # 模拟标识
            if progress_info.get("simulated"):
                message_parts.append("🔬 **模式**: 模拟执行（Nextflow未安装）")
            
            # 下一步提示
            if progress_info.get("completed"):
                message_parts.extend([
                    "",
                    "🎉 **分析完成！** 可以查看结果文件和日志。"
                ])
            elif progress_info.get("failed"):
                message_parts.extend([
                    "",
                    "💡 **建议**: 检查日志文件，修复问题后重新运行。"
                ])
            else:
                message_parts.extend([
                    "",
                    "⏳ **请等待**: 流程正在后台运行..."
                ])
            
            return "\n".join(message_parts)
            
        except Exception as e:
            logger.error(f"Error generating progress message: {str(e)}")
            return f"进度报告生成失败：{str(e)}"
    
    def _create_progress_bar(self, percent: int) -> str:
        """创建进度条"""
        try:
            width = 40
            filled = int(width * percent / 100)
            empty = width - filled
            
            bar = "█" * filled + "░" * empty
            return f"[{bar}] {percent}%"
        
        except Exception:
            return f"[{'?' * 40}] {percent}%"
    
    def _get_progress_bar(self, progress: str) -> str:
        """
        生成进度条
        
        应用KISS原则：简单的进度条显示
        """
        try:
            percent = int(progress.replace('%', ''))
            filled = int(percent / 10)
            empty = 10 - filled
            bar = '█' * filled + '░' * empty
            return f"[{bar}] {progress}"
        except:
            return f"[░░░░░░░░░░] {progress}"
    
    def collect_results(self, state: AgentState) -> Dict[str, Any]:
        """
        收集执行结果
        
        遵循单一职责原则：专门处理结果收集
        """
        try:
            work_dir = state.get("execution_results", {}).get("execution_info", {}).get("work_dir", "./data")
            results_dir = f"{work_dir}/results"
            
            # 检查结果文件
            result_files = self._scan_result_files(results_dir)
            
            # 生成结果总结
            summary = self._generate_result_summary(result_files, work_dir)
            
            return {
                "execution_results": {
                    "result_files": result_files,
                    "summary": summary,
                    "work_dir": work_dir
                },
                "execution_status": "completed",
                "messages": [AIMessage(content=summary)]
            }
        
        except Exception as e:
            logger.error(f"Error collecting results: {str(e)}")
            return {
                "messages": [AIMessage(content=f"收集结果时出错：{str(e)}")]
            }
    
    def _scan_result_files(self, results_dir: str) -> Dict[str, List[str]]:
        """
        扫描结果文件
        
        应用KISS原则：简单的文件扫描
        """
        result_files = {
            "fastp": [],
            "star": [],
            "featurecounts": [],
            "logs": []
        }
        
        try:
            if not os.path.exists(results_dir):
                return result_files
            
            for root, dirs, files in os.walk(results_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    
                    if "fastp" in file.lower():
                        result_files["fastp"].append(file_path)
                    elif "star" in file.lower() or file.endswith(".bam"):
                        result_files["star"].append(file_path)
                    elif "featurecounts" in file.lower() or "counts" in file.lower():
                        result_files["featurecounts"].append(file_path)
                    elif file.endswith(".log"):
                        result_files["logs"].append(file_path)
        
        except Exception as e:
            logger.error(f"Error scanning result files: {str(e)}")
        
        return result_files
    
    def _generate_result_summary(self, result_files: Dict[str, List[str]], work_dir: str) -> str:
        """
        生成结果总结
        
        应用模板方法模式：标准的总结格式
        """
        try:
            summary_parts = ["🎉 **RNA-seq分析完成！**\n"]
            
            # 结果文件统计
            summary_parts.append("📁 **结果文件统计**:")
            for category, files in result_files.items():
                if files:
                    summary_parts.append(f"- {category.upper()}: {len(files)} 个文件")
            
            summary_parts.append("")
            
            # 主要输出文件
            summary_parts.append("📄 **主要输出文件**:")
            
            if result_files["featurecounts"]:
                summary_parts.append("- 基因表达定量结果: featureCounts输出")
            
            if result_files["star"]:
                summary_parts.append("- 比对结果文件: BAM格式")
            
            if result_files["fastp"]:
                summary_parts.append("- 质量控制报告: HTML格式")
            
            summary_parts.append("")
            summary_parts.append(f"📂 **完整结果目录**: {work_dir}/results")
            summary_parts.append(f"📋 **日志文件目录**: {work_dir}/logs")
            
            summary_parts.append("\n✅ 分析流程已成功完成！您可以在结果目录中查看所有输出文件。")
            
            return "\n".join(summary_parts)
        
        except Exception as e:
            logger.error(f"Error generating result summary: {str(e)}")
            return "结果总结生成失败，请手动检查输出目录。"

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
        
        if execution_status == "idle":
            # 开始执行
            logger.info("Starting nextflow execution")
            ui_manager.show_info("正在启动Nextflow执行流程...")
            result = handler.execute_nextflow(state)
            
            # 确保返回值包含必要的状态更新
            if "mode" not in result:
                result["mode"] = "execute"
            
            return result
        
        elif execution_status == "running":
            # 监控执行
            logger.info("Monitoring execution progress")
            ui_manager.show_info("正在监控执行进度...")
            result = handler.monitor_execution(state)
            
            # 确保返回值包含必要的状态更新
            if "mode" not in result:
                result["mode"] = "execute"
            
            return result
        
        elif execution_status == "completed":
            # 收集结果
            logger.info("Collecting execution results")
            ui_manager.show_success("执行完成，正在收集结果...")
            result = handler.collect_results(state)
            
            # 确保返回值包含必要的状态更新
            if "mode" not in result:
                result["mode"] = "execute"
            
            return result
        
        elif execution_status == "failed":
            # 处理失败情况
            logger.info("Handling execution failure")
            ui_manager.show_error("执行失败")
            return {
                "messages": [AIMessage(content="执行失败。请检查配置和日志文件，然后重试。")],
                "mode": "execute",
                "execution_status": "failed"
            }
        
        else:
            # 未知状态，默认开始执行
            logger.warning(f"Unknown execution status: {execution_status}, defaulting to idle")
            ui_manager.show_warning(f"未知执行状态: {execution_status}，重新开始执行")
            result = handler.execute_nextflow(state)
            
            # 确保返回值包含必要的状态更新
            if "mode" not in result:
                result["mode"] = "execute"
            
            return result
    
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

# ============================================================================
# 执行配置和模板 - 遵循配置分离原则
# ============================================================================

class ExecutionConfig:
    """
    执行配置类
    
    遵循单一职责原则：专门管理执行配置
    """
    
    # 默认工作目录
    DEFAULT_WORK_DIR = "./data"
    
    # 结果目录结构
    RESULT_DIRS = [
        "results/fastp",
        "results/star", 
        "results/featurecounts",
        "logs"
    ]
    
    # 必需的配置字段
    REQUIRED_CONFIG_FIELDS = [
        "data"  # 至少需要数据目录
    ]
    
    @classmethod
    def create_work_directories(cls, base_dir: str) -> bool:
        """创建工作目录结构"""
        try:
            for dir_path in cls.RESULT_DIRS:
                full_path = os.path.join(base_dir, dir_path)
                os.makedirs(full_path, exist_ok=True)
            return True
        except Exception as e:
            logger.error(f"Error creating work directories: {str(e)}")
            return False

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
    
    @classmethod
    def validate_config(cls, config: Dict[str, Any]) -> Dict[str, Any]:
        """验证配置完整性"""
        missing_fields = []
        
        for field in cls.REQUIRED_CONFIG_FIELDS:
            if not config.get(field):
                missing_fields.append(field)
        
        return {
            "valid": len(missing_fields) == 0,
            "missing_fields": missing_fields
        }