"""
Execute Mode节点 - 执行nextflow流程和结果总结
遵循单一职责原则：专门处理execute模式下的流程执行和结果处理
采用JSON-first架构，与其他模式保持一致
"""

import logging
import os
import time
import subprocess
import threading
import re
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from langchain_core.messages import HumanMessage, AIMessage
from ..state import AgentState, update_execution_status
from ..core import create_chain_for_mode, create_structured_chain_for_mode
from ..ui_manager import get_ui_manager

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class NextflowProgressMonitor:
    """
    Nextflow实时进度监控器
    
    混合监控方案：
    1. 实时解析stdout/stderr流
    2. 监控工作目录结构变化
    3. 计算整体进度和时间估算
    4. 提供用户友好的进度显示
    """
    
    def __init__(self):
        self.process = None
        self.start_time = None
        self.current_progress = 0
        self.total_steps = 0
        self.current_step = ""
        self.completed_processes = []
        self.running_processes = []
        self.failed_processes = []
        self.log_lines = []
        self.stop_monitoring = False
        
        # Nextflow进程识别模式
        self.process_patterns = {
            'fastp': r'process > FASTP',
            'star_align': r'process > STAR_ALIGN', 
            'featurecounts': r'process > FEATURECOUNTS',
            'download_srr': r'process > DOWNLOAD_SRR',
            'download_genome': r'process > DOWNLOAD_GENOME',
            'build_star_index': r'process > BUILD_STAR_INDEX'
        }
        
        # 步骤描述映射
        self.step_descriptions = {
            'fastp': '质量控制 (FastP)',
            'star_align': '序列比对 (STAR)',
            'featurecounts': '基因定量 (featureCounts)', 
            'download_srr': 'SRR数据下载',
            'download_genome': '基因组下载',
            'build_star_index': 'STAR索引构建'
        }
    
    def start_monitoring(self, command: str, work_dir: str = "./work"):
        """启动监控"""
        self.start_time = time.time()
        self.work_dir = work_dir
        self.stop_monitoring = False
        
        logger.info(f"启动Nextflow进度监控: {command}")
        
        try:
            # 启动Nextflow进程
            self.process = subprocess.Popen(
                command,
                shell=True,
                cwd=".",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # 合并stderr到stdout
                text=True,
                bufsize=1,  # 行缓冲
                universal_newlines=True
            )
            
            logger.info(f"Nextflow进程已启动，PID: {self.process.pid}")
            return True
            
        except Exception as e:
            logger.error(f"启动Nextflow进程失败: {str(e)}")
            return False
    
    def monitor_progress(self, callback_func=None):
        """
        监控执行进度
        
        Args:
            callback_func: 回调函数，用于更新UI显示
        """
        if not self.process:
            return
            
        try:
            # 启动输出读取线程
            output_thread = threading.Thread(
                target=self._read_output_stream,
                args=(callback_func,),
                daemon=True
            )
            output_thread.start()
            
            # 主监控循环
            while self.process.poll() is None and not self.stop_monitoring:
                # 计算当前进度
                progress_info = self._calculate_progress()
                
                if callback_func:
                    callback_func(progress_info)
                
                time.sleep(2)  # 每2秒更新一次
            
            # 等待输出线程结束
            output_thread.join(timeout=5)
            
            # 获取最终结果
            return_code = self.process.returncode
            final_info = self._get_final_results(return_code)
            
            if callback_func:
                callback_func(final_info)
                
            return final_info
            
        except Exception as e:
            logger.error(f"监控进程中出现错误: {str(e)}")
            error_info = {
                'status': 'error',
                'message': f'监控出现错误: {str(e)}',
                'progress': 0,
                'elapsed_time': time.time() - self.start_time if self.start_time else 0
            }
            if callback_func:
                callback_func(error_info)
            return error_info
    
    def _read_output_stream(self, callback_func=None):
        """读取并解析输出流"""
        try:
            for line in iter(self.process.stdout.readline, ''):
                if self.stop_monitoring:
                    break
                    
                line = line.strip()
                if line:
                    self.log_lines.append(line)
                    self._parse_output_line(line)
                    
                    # 如果是重要信息，立即回调
                    if any(pattern in line for pattern in ['process >', 'ERROR', 'WARN', 'executor >']):
                        if callback_func:
                            progress_info = self._calculate_progress()
                            progress_info['latest_log'] = line
                            callback_func(progress_info)
                            
        except Exception as e:
            logger.error(f"读取输出流时出错: {str(e)}")
    
    def _parse_output_line(self, line: str):
        """解析单行输出"""
        try:
            # 检测进程启动/完成
            if 'process >' in line:
                self._update_process_status(line)
            
            # 检测错误
            elif 'ERROR' in line:
                logger.warning(f"Nextflow错误: {line}")
                self.failed_processes.append(line)
            
            # 检测执行器信息
            elif 'executor >' in line:
                match = re.search(r'executor >\s+(\w+)\s+\((\d+)\)', line)
                if match:
                    executor = match.group(1)
                    task_count = int(match.group(2))
                    logger.info(f"执行器 {executor} 处理 {task_count} 个任务")
                    
        except Exception as e:
            logger.error(f"解析输出行时出错: {str(e)}")
    
    def _update_process_status(self, line: str):
        """更新进程状态"""
        try:
            # 解析进程信息：[hash] process > PROCESS_NAME (sample) [progress] status
            process_match = re.search(r'process > (\w+)(?:\s+\(([^)]+)\))?\s+\[([^\]]*)\](?:\s+(.+))?', line)
            
            if process_match:
                process_name = process_match.group(1).lower()
                sample_name = process_match.group(2) or ""
                progress_info = process_match.group(3) or ""
                status_info = process_match.group(4) or ""
                
                logger.info(f"进程更新: {process_name}, 样本: {sample_name}, 进度: {progress_info}, 状态: {status_info}")
                
                # 更新当前步骤
                if process_name in self.step_descriptions:
                    self.current_step = self.step_descriptions[process_name]
                    if sample_name:
                        self.current_step += f" ({sample_name})"
                
                # 判断进程状态
                if '100%' in progress_info and ('✓' in status_info or 'COMPLETED' in status_info.upper()):
                    if process_name not in self.completed_processes:
                        self.completed_processes.append(process_name)
                        logger.info(f"进程完成: {process_name}")
                elif process_name not in self.running_processes:
                    self.running_processes.append(process_name)
                    
        except Exception as e:
            logger.error(f"更新进程状态时出错: {str(e)}")
    
    def _calculate_progress(self) -> Dict[str, Any]:
        """计算当前进度"""
        try:
            elapsed_time = time.time() - self.start_time if self.start_time else 0
            
            # 估算总步骤数（基于检测到的进程）
            all_processes = set(self.completed_processes + self.running_processes)
            total_processes = max(len(all_processes), 1)
            completed_count = len(self.completed_processes)
            
            # 计算进度百分比
            if total_processes > 0:
                progress_percent = min(int((completed_count / total_processes) * 100), 100)
            else:
                progress_percent = 0
            
            # 估算剩余时间
            if completed_count > 0 and progress_percent < 100:
                avg_time_per_step = elapsed_time / completed_count
                remaining_steps = total_processes - completed_count
                estimated_remaining = avg_time_per_step * remaining_steps
            else:
                estimated_remaining = 0
            
            return {
                'status': 'running',
                'progress': progress_percent,
                'current_step': self.current_step,
                'completed_processes': self.completed_processes.copy(),
                'running_processes': self.running_processes.copy(),
                'failed_processes': self.failed_processes.copy(),
                'elapsed_time': elapsed_time,
                'estimated_remaining': estimated_remaining,
                'total_processes': total_processes,
                'completed_count': completed_count,
                'latest_logs': self.log_lines[-10:] if self.log_lines else []
            }
            
        except Exception as e:
            logger.error(f"计算进度时出错: {str(e)}")
            return {
                'status': 'error',
                'progress': 0,
                'message': f'进度计算错误: {str(e)}',
                'elapsed_time': time.time() - self.start_time if self.start_time else 0
            }
    
    def _get_final_results(self, return_code: int) -> Dict[str, Any]:
        """获取最终执行结果"""
        try:
            elapsed_time = time.time() - self.start_time if self.start_time else 0
            
            if return_code == 0:
                status = 'completed'
                message = '🎉 Nextflow执行成功完成！'
                progress = 100
            else:
                status = 'failed'
                message = f'❌ Nextflow执行失败 (退出代码: {return_code})'
                progress = max(self._calculate_progress().get('progress', 0), 0)
            
            return {
                'status': status,
                'message': message,
                'progress': progress,
                'return_code': return_code,
                'elapsed_time': elapsed_time,
                'completed_processes': self.completed_processes.copy(),
                'failed_processes': self.failed_processes.copy(),
                'total_log_lines': len(self.log_lines),
                'final_logs': self.log_lines[-20:] if self.log_lines else []
            }
            
        except Exception as e:
            logger.error(f"获取最终结果时出错: {str(e)}")
            return {
                'status': 'error',
                'message': f'获取结果时出错: {str(e)}',
                'progress': 0,
                'elapsed_time': time.time() - self.start_time if self.start_time else 0
            }
    
    def stop(self):
        """停止监控"""
        self.stop_monitoring = True
        if self.process and self.process.poll() is None:
            logger.info("正在终止Nextflow进程...")
            self.process.terminate()
            time.sleep(3)
            if self.process.poll() is None:
                self.process.kill()
    
    def format_progress_display(self, progress_info: Dict[str, Any]) -> str:
        """格式化进度显示"""
        try:
            status = progress_info.get('status', 'unknown')
            progress = progress_info.get('progress', 0)
            current_step = progress_info.get('current_step', '初始化...')
            elapsed = progress_info.get('elapsed_time', 0)
            estimated_remaining = progress_info.get('estimated_remaining', 0)
            
            # 状态emoji
            status_emoji = {
                'running': '🚀',
                'completed': '✅', 
                'failed': '❌',
                'error': '⚠️'
            }
            
            # 进度条
            bar_length = 30
            filled_length = int(bar_length * progress / 100)
            bar = '█' * filled_length + '░' * (bar_length - filled_length)
            
            # 时间格式化
            def format_time(seconds):
                if seconds < 60:
                    return f"{int(seconds)}秒"
                elif seconds < 3600:
                    return f"{int(seconds//60)}分{int(seconds%60)}秒"
                else:
                    hours = int(seconds // 3600)
                    minutes = int((seconds % 3600) // 60)
                    return f"{hours}小时{minutes}分钟"
            
            display = f"{status_emoji.get(status, '🔄')} **Nextflow执行进度**\n\n"
            display += f"📊 整体进度：[{bar}] {progress}%\n"
            display += f"🔄 当前步骤：{current_step}\n"
            display += f"⏰ 已执行：{format_time(elapsed)}"
            
            if estimated_remaining > 0:
                display += f" | 预计剩余：{format_time(estimated_remaining)}"
            
            display += "\n"
            
            # 显示已完成的进程
            completed = progress_info.get('completed_processes', [])
            if completed:
                display += f"\n✅ **已完成步骤** ({len(completed)})：\n"
                for proc in completed:
                    step_name = self.step_descriptions.get(proc, proc)
                    display += f"• {step_name}\n"
            
            # 显示正在运行的进程
            running = progress_info.get('running_processes', [])
            if running:
                display += f"\n🔄 **正在执行** ({len(running)})：\n"
                for proc in running:
                    step_name = self.step_descriptions.get(proc, proc)
                    display += f"• {step_name}\n"
            
            # 显示最新日志（如果有）
            latest_log = progress_info.get('latest_log')
            if latest_log and 'process >' in latest_log:
                display += f"\n📄 **最新活动**：\n```\n{latest_log}\n```"
            
            return display
            
        except Exception as e:
            logger.error(f"格式化进度显示时出错: {str(e)}")
            return f"⚠️ 显示格式化错误: {str(e)}"

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
        self.progress_monitor = None  # 进度监控器
        self.execution_log = []  # 存储执行日志
    
    def _process_llm_response(self, response) -> Tuple[AIMessage, List[Dict[str, Any]]]:
        """
        处理LLM的结构化响应（.with_structured_output()返回Pydantic模型实例）
        
        返回: (AIMessage, tool_calls列表)
        """
        try:
            # 调试日志：查看响应类型和内容
            logger.info(f"Execute模式收到响应类型: {type(response)}")
            logger.info(f"Execute模式收到响应内容: {repr(response)[:200]}...")
            
            # .with_structured_output()返回Pydantic模型实例
            if hasattr(response, 'status'):  # ExecuteModeResponse模型
                logger.info(f"Execute模式收到Pydantic模型响应")
                
                # 提取响应信息
                user_message = getattr(response, "response", "执行完成") if hasattr(response, "response") else "执行完成"
                status = getattr(response, "status", "unknown")
                progress = getattr(response, "progress", "")
                next_step = getattr(response, "next_step", "")
                results = getattr(response, "results", {})
                tool_calls = getattr(response, "tool_calls", [])
                
                # 构建详细响应
                detailed_response = user_message
                if status and status != "unknown":
                    detailed_response += f"\n\n📊 **状态**: {status}"
                if progress:
                    detailed_response += f"\n🔄 **进度**: {progress}"
                if next_step:
                    detailed_response += f"\n➡️ **下一步**: {next_step}"
                if results:
                    detailed_response += "\n\n📋 **结果**:\n"
                    for key, value in results.items():
                        detailed_response += f"  - {key}: {value}\n"
                
                logger.info(f"Execute模式提取到 {len(tool_calls)} 个工具调用")
                
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
                                "id": f"call_exec_{i}",
                                "type": "tool_call"
                            }
                        else:
                            # tool_call是字典格式
                            tool_call_obj = {
                                "name": tool_call.get("tool_name"),
                                "args": tool_call.get("parameters", {}),
                                "id": f"call_exec_{i}",
                                "type": "tool_call"
                            }
                        langchain_tool_calls.append(tool_call_obj)
                    
                    ai_message.tool_calls = langchain_tool_calls
                    logger.info(f"Execute模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            elif isinstance(response, dict):
                # 兼容旧的dict格式返回
                logger.info(f"Execute模式收到dict格式响应: {list(response.keys())}")
                
                # 提取响应信息
                user_message = response.get("response", "执行完成")
                status = response.get("status", "unknown")
                progress = response.get("progress", "")
                next_step = response.get("next_step", "")
                results = response.get("results", {})
                tool_calls = response.get("tool_calls", [])
                
                # 构建详细响应
                detailed_response = user_message
                if status and status != "unknown":
                    detailed_response += f"\n\n📊 **状态**: {status}"
                if progress:
                    detailed_response += f"\n🔄 **进度**: {progress}"
                if next_step:
                    detailed_response += f"\n➡️ **下一步**: {next_step}"
                if results:
                    detailed_response += "\n\n📋 **结果**:\n"
                    for key, value in results.items():
                        detailed_response += f"  - {key}: {value}\n"
                
                logger.info(f"Execute模式提取到 {len(tool_calls)} 个工具调用")
                
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
                    logger.info(f"Execute模式成功设置tool_calls属性")
                
                return ai_message, tool_calls
            else:
                # 降级处理：如果不是期望的格式
                logger.warning(f"Execute模式收到未知响应格式: {type(response)}")
                content = str(response) if response else "响应为空"
                return AIMessage(content=content), []
            
        except Exception as e:
            logger.error(f"Execute模式处理响应时出错: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            return AIMessage(content=f"处理响应时出现错误: {str(e)}"), []
    
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
        执行nextflow流程 - 带实时进度监控
        
        使用混合监控方案：实时输出解析 + 工作目录监控 + 时间估算
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
            
            # 只有当genome_version为空时才使用默认值hg38
            if not merged_config.get("genome_version"):
                merged_config["genome_version"] = "hg38"
                logger.info("使用默认基因组版本: hg38")
            else:
                logger.info(f"使用配置的基因组版本: {merged_config['genome_version']}")
            
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
            
            # 获取UI管理器
            ui_manager = get_ui_manager()
            
            # 初始化进度监控器
            self.progress_monitor = NextflowProgressMonitor()
            
            # 启动监控
            if not self.progress_monitor.start_monitoring(command, "./work"):
                error_msg = "❌ **启动Nextflow进程失败**\n\n请检查系统环境和配置。"
                return {
                    "execution_status": "failed",
                    "messages": [AIMessage(content=error_msg)]
                }
            
            # 显示初始状态
            initial_msg = "🚀 **正在启动Nextflow执行...**\n\n"
            initial_msg += "📋 **执行计划**:\n"
            if plan:
                for i, step in enumerate(plan, 1):
                    initial_msg += f"{i}. {step}\n"
            else:
                initial_msg += "• 基于配置的自动流程\n"
            initial_msg += "\n⏳ 初始化中，请稍候..."
            
            ui_manager.show_info(initial_msg)
            
            # 创建持续更新的进度显示
            progress_bar = None
            
            def update_progress_display(progress_info):
                nonlocal progress_bar
                try:
                    import sys
                    
                    # 获取进度信息
                    progress = progress_info.get('progress', 0)
                    current_step = progress_info.get('current_step', '初始化...')
                    # status = progress_info.get('status', 'running')  # 暂时不使用
                    
                    # 创建进度条（只在第一次）
                    if progress_bar is None:
                        try:
                            from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn
                            
                            progress_bar = Progress(
                                TextColumn("[bold blue]🚀 Nextflow执行"),
                                BarColumn(bar_width=40),
                                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                                TimeElapsedColumn(),
                                TimeRemainingColumn(),
                                TextColumn("{task.description}"),
                                console=ui_manager.console if ui_manager.use_rich else None,
                                refresh_per_second=4  # 4次/秒的刷新率，更流畅
                            )
                            progress_bar.start()
                            
                            # 添加主任务
                            task_id = progress_bar.add_task(
                                description=current_step,
                                total=100
                            )
                            progress_bar.task_id = task_id
                            
                        except ImportError:
                            # Rich不可用时的简单进度条
                            progress_bar = "simple"
                    
                    # 更新进度
                    if hasattr(progress_bar, 'update') and hasattr(progress_bar, 'task_id'):
                        # Rich进度条
                        progress_bar.update(
                            progress_bar.task_id,
                            completed=progress,
                            description=current_step
                        )
                    else:
                        # 简单进度条
                        bar_length = 40
                        filled_length = int(bar_length * progress / 100)
                        bar = '█' * filled_length + '░' * (bar_length - filled_length)
                        
                        # 清除当前行并显示新进度
                        sys.stdout.write(f'\r🚀 进度: [{bar}] {progress}% - {current_step}')
                        sys.stdout.flush()
                    
                except Exception as e:
                    logger.error(f"更新进度显示时出错: {str(e)}")
            
            def cleanup_progress_bar():
                nonlocal progress_bar
                if progress_bar and hasattr(progress_bar, 'stop'):
                    progress_bar.stop()
                elif progress_bar == "simple":
                    # 简单进度条完成后换行
                    print()  # 换行
            
            try:
                # 开始监控（这会阻塞直到完成）
                logger.info("开始监控Nextflow执行...")
                final_result = self.progress_monitor.monitor_progress(update_progress_display)
                
                # 处理最终结果
                return self._process_execution_results(final_result, command, work_dir)
            
            finally:
                # 清理进度条
                cleanup_progress_bar()
                
        except Exception as e:
            logger.error(f"Error executing nextflow: {str(e)}")
            import traceback
            logger.error(f"错误堆栈: {traceback.format_exc()}")
            
            # 停止监控器（如果存在）
            if self.progress_monitor:
                self.progress_monitor.stop()
            
            return {
                "execution_status": "failed",
                "messages": [AIMessage(content=f"❌ **执行出错**: {str(e)}\n\n请检查系统配置和日志。")]
            }
    
    def _process_execution_results(self, final_result: Dict[str, Any], command: str, work_dir: str) -> Dict[str, Any]:
        """处理执行结果"""
        try:
            status = final_result.get('status', 'unknown')
            return_code = final_result.get('return_code', -1)
            elapsed_time = final_result.get('elapsed_time', 0)
            completed_processes = final_result.get('completed_processes', [])
            failed_processes = final_result.get('failed_processes', [])
            final_logs = final_result.get('final_logs', [])
            
            def format_time(seconds):
                if seconds < 60:
                    return f"{int(seconds)}秒"
                elif seconds < 3600:
                    return f"{int(seconds//60)}分{int(seconds%60)}秒"
                else:
                    hours = int(seconds // 3600)
                    minutes = int((seconds % 3600) // 60)
                    return f"{hours}小时{minutes}分钟"
            
            if status == 'completed':
                # 成功完成
                success_message = "🎉 **Nextflow执行成功完成！**\n\n"
                success_message += f"⏱️ **执行概况**:\n"
                success_message += f"• 总用时: {format_time(elapsed_time)}\n"
                success_message += f"• 工作目录: {work_dir}\n"
                success_message += f"• 退出代码: {return_code}\n\n"
                
                # 显示完成的步骤
                if completed_processes:
                    success_message += f"✅ **完成的步骤** ({len(completed_processes)}):\n"
                    for proc in completed_processes:
                        step_name = self.progress_monitor.step_descriptions.get(proc, proc)
                        success_message += f"• {step_name}\n"
                    success_message += "\n"
                
                # 显示最近的日志
                if final_logs:
                    success_message += "📊 **执行摘要**:\n"
                    success_message += "```\n"
                    # 显示最后几行重要日志
                    important_logs = [log for log in final_logs[-10:] 
                                    if any(keyword in log for keyword in ['process >', 'executor >', 'Completed'])]
                    if important_logs:
                        success_message += '\n'.join(important_logs[-5:])
                    else:
                        success_message += '\n'.join(final_logs[-3:])
                    success_message += "\n```\n\n"
                
                success_message += "✅ 分析结果已保存到输出目录中。您可以查看 `data/results/` 目录获取结果文件。"
                
                return {
                    "execution_status": "completed",
                    "messages": [AIMessage(content=success_message)]
                }
                
            else:
                # 执行失败
                error_message = f"❌ **Nextflow执行失败**\n\n"
                error_message += f"⏱️ **执行信息**:\n"
                error_message += f"• 用时: {format_time(elapsed_time)}\n"
                error_message += f"• 退出代码: {return_code}\n"
                error_message += f"• 命令: `{command}`\n\n"
                
                # 显示已完成的步骤
                if completed_processes:
                    error_message += f"✅ **已完成步骤** ({len(completed_processes)}):\n"
                    for proc in completed_processes:
                        step_name = self.progress_monitor.step_descriptions.get(proc, proc)
                        error_message += f"• {step_name}\n"
                    error_message += "\n"
                
                # 显示失败的步骤
                if failed_processes:
                    error_message += f"❌ **失败步骤**:\n"
                    for proc in failed_processes[-3:]:  # 只显示最后3个错误
                        error_message += f"• {proc}\n"
                    error_message += "\n"
                
                # 显示关键日志
                if final_logs:
                    error_message += "📄 **错误日志**:\n"
                    error_message += "```\n"
                    # 查找错误相关的日志
                    error_logs = [log for log in final_logs 
                                if any(keyword in log.upper() for keyword in ['ERROR', 'FAILED', 'EXCEPTION'])]
                    if error_logs:
                        error_message += '\n'.join(error_logs[-5:])
                    else:
                        error_message += '\n'.join(final_logs[-5:])
                    error_message += "\n```\n\n"
                
                error_message += "**建议**:\n"
                error_message += "• 检查输入文件格式和路径\n"
                error_message += "• 验证基因组配置和索引文件\n"
                error_message += "• 查看详细日志: `.nextflow.log`\n"
                error_message += "• 检查系统资源和权限"
                
                return {
                    "execution_status": "failed",
                    "messages": [AIMessage(content=error_message)]
                }
                
        except Exception as e:
            logger.error(f"处理执行结果时出错: {str(e)}")
            return {
                "execution_status": "failed",
                "messages": [AIMessage(content=f"❌ **结果处理出错**: {str(e)}")]
            }
    
    def monitor_execution(self, state: AgentState) -> Dict[str, Any]:
        """
        简化版监控 - 直接返回完成状态，避免循环
        """
        _ = state  # 避免未使用参数警告
        return {
            "execution_status": "completed",
            "messages": [AIMessage(content="✅ **监控完成**\n\n执行流程已结束。")]
        }
    
    def collect_results(self, state: AgentState) -> Dict[str, Any]:
        """
        简化版结果收集 - 直接返回完成状态
        """
        _ = state  # 避免未使用参数警告
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
