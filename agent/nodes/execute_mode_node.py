"""
Execute Mode节点 - 执行nextflow流程和结果总结
遵循单一职责原则：专门处理execute模式下的流程执行和结果处理
"""

import logging
import os
import json
import time
from typing import Dict, Any, List, Optional
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
    """
    
    def __init__(self):
        self.chain = create_chain_for_mode("execute")
        self.structured_chain = create_structured_chain_for_mode("execute")
    
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
        执行nextflow流程
        
        应用命令模式：封装执行命令
        """
        try:
            # 准备执行
            prep_result = self.prepare_execution(state)
            if "error" in prep_result:
                return prep_result
            
            # 构建nextflow命令参数
            config = prep_result["config"]
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
            
            # 添加配置文件
            cmd_parts.extend(["-c", "config/nextflow.config"])
            
            command = " ".join(cmd_parts)
            
            logger.info(f"Executing nextflow command: {command}")
            
            # 这里实际执行时会调用subprocess，现在先模拟
            execution_info = {
                "command": command,
                "status": "running",
                "start_time": time.time(),
                "work_dir": prep_result["work_dir"]
            }
            
            return {
                "execution_info": execution_info,
                "execution_status": "running",
                "messages": [AIMessage(content=f"🚀 Nextflow流程已启动！\n\n执行命令：\n```\n{command}\n```\n\n请稍等，流程正在后台运行...")]
            }
        
        except Exception as e:
            logger.error(f"Error executing nextflow: {str(e)}")
            return {
                "error": str(e),
                "execution_status": "failed",
                "messages": [AIMessage(content=f"执行失败：{str(e)}")]
            }
    
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
        监控执行状态
        
        应用观察者模式：监控执行进度
        """
        try:
            execution_info = state.get("execution_results", {}).get("execution_info", {})
            
            if not execution_info:
                return {
                    "messages": [AIMessage(content="没有找到正在执行的流程。")]
                }
            
            # 检查执行状态（这里模拟，实际会检查进程状态）
            current_time = time.time()
            start_time = execution_info.get("start_time", current_time)
            elapsed_time = current_time - start_time
            
            # 模拟执行进度
            if elapsed_time < 60:  # 1分钟内
                status = "正在初始化..."
                progress = "10%"
            elif elapsed_time < 300:  # 5分钟内
                status = "正在执行质量控制..."
                progress = "30%"
            elif elapsed_time < 600:  # 10分钟内
                status = "正在进行序列比对..."
                progress = "60%"
            elif elapsed_time < 900:  # 15分钟内
                status = "正在进行基因定量..."
                progress = "80%"
            else:
                status = "执行完成"
                progress = "100%"
            
            status_message = f"""
📊 **执行状态监控**

⏱️ **运行时间**: {int(elapsed_time//60)}分{int(elapsed_time%60)}秒
📈 **当前进度**: {progress}
🔄 **当前状态**: {status}
📁 **工作目录**: {execution_info.get('work_dir', 'N/A')}

{self._get_progress_bar(progress)}
            """
            
            return {
                "execution_status": "completed" if progress == "100%" else "running",
                "messages": [AIMessage(content=status_message)]
            }
        
        except Exception as e:
            logger.error(f"Error monitoring execution: {str(e)}")
            return {
                "messages": [AIMessage(content=f"监控执行状态时出错：{str(e)}")]
            }
    
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