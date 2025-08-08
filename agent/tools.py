import os
import json
import glob
import subprocess
from typing import Dict, Any, List, Optional
from pathlib import Path
from langchain_core.tools import tool
from pydantic.v1 import BaseModel, Field

# ============================================================================
# 输入模型定义 - 遵循接口隔离原则
# ============================================================================

class DirectoryQueryArgs(BaseModel):
    """目录查询参数模型"""
    directory_path: str = Field(description="要查询的目录路径")

class FastqQueryArgs(BaseModel):
    """FASTQ文件查询参数模型"""
    directory_path: str = Field(description="包含FASTQ文件的目录路径")
    pattern: str = Field(default="*.fastq*", description="文件匹配模式")

class GenomeQueryArgs(BaseModel):
    """基因组查询参数模型"""
    genome_name: str = Field(description="基因组名称，如hg38、mm39等")
    config_path: str = Field(default="config/genomes.json", description="基因组配置文件路径")

class NextflowConfigArgs(BaseModel):
    """Nextflow配置参数模型"""
    param_name: str = Field(description="参数名称")
    param_value: Any = Field(description="参数值")

class BatchConfigArgs(BaseModel):
    """批量配置更新参数模型"""
    config_updates: Dict[str, Any] = Field(description="要更新的配置字典")

class ModeSwitch(BaseModel):
    """模式切换参数模型"""
    target_mode: str = Field(description="目标模式：plan或execute")
    reason: str = Field(description="切换原因")

class ExecutionArgs(BaseModel):
    """执行参数模型"""
    config_path: str = Field(default="config/nextflow.config", description="nextflow配置文件路径")
    work_dir: str = Field(default="./work", description="工作目录")

# ============================================================================
# 信息查询工具组 - 遵循单一职责原则
# ============================================================================

@tool(args_schema=DirectoryQueryArgs)
def list_directory_contents(directory_path: str) -> str:
    """
    列出指定目录的内容
    
    应用KISS原则：简单直接的目录列表功能
    """
    try:
        if not os.path.exists(directory_path):
            return f"错误：目录 '{directory_path}' 不存在"
        
        if not os.path.isdir(directory_path):
            return f"错误：'{directory_path}' 不是一个目录"
        
        contents = []
        for item in os.listdir(directory_path):
            item_path = os.path.join(directory_path, item)
            if os.path.isdir(item_path):
                contents.append(f"📁 {item}/")
            else:
                size = os.path.getsize(item_path)
                contents.append(f"📄 {item} ({size} bytes)")
        
        if not contents:
            return f"目录 '{directory_path}' 为空"
        
        return f"目录 '{directory_path}' 内容：\n" + "\n".join(contents)
    
    except PermissionError:
        return f"错误：没有权限访问目录 '{directory_path}'"
    except Exception as e:
        return f"查询目录时发生错误：{str(e)}"

@tool(args_schema=FastqQueryArgs)
def query_fastq_files(directory_path: str, pattern: str = "*.fastq*") -> str:
    """
    查询指定目录下的FASTQ文件信息
    
    遵循DRY原则：统一的FASTQ文件查询逻辑
    """
    try:
        if not os.path.exists(directory_path):
            return f"错误：目录 '{directory_path}' 不存在"
        
        # 使用glob查找FASTQ文件
        search_pattern = os.path.join(directory_path, "**", pattern)
        fastq_files = glob.glob(search_pattern, recursive=True)
        
        if not fastq_files:
            return f"在目录 '{directory_path}' 中未找到匹配 '{pattern}' 的FASTQ文件"
        
        # 分析文件信息
        paired_files = {}
        single_files = []
        
        for file_path in sorted(fastq_files):
            file_name = os.path.basename(file_path)
            file_size = os.path.getsize(file_path)
            
            # 判断是否为双端测序文件
            if "_1.fastq" in file_name or "_R1" in file_name:
                sample_id = file_name.replace("_1.fastq", "").replace("_R1", "").split(".")[0]
                if sample_id not in paired_files:
                    paired_files[sample_id] = {}
                paired_files[sample_id]["R1"] = {"path": file_path, "size": file_size}
            elif "_2.fastq" in file_name or "_R2" in file_name:
                sample_id = file_name.replace("_2.fastq", "").replace("_R2", "").split(".")[0]
                if sample_id not in paired_files:
                    paired_files[sample_id] = {}
                paired_files[sample_id]["R2"] = {"path": file_path, "size": file_size}
            else:
                single_files.append({"name": file_name, "path": file_path, "size": file_size})
        
        # 构建结果
        result = [f"FASTQ文件查询结果 (目录: {directory_path})：\n"]
        
        if paired_files:
            result.append("双端测序文件：")
            for sample_id, files in paired_files.items():
                result.append(f"  样本: {sample_id}")
                if "R1" in files:
                    result.append(f"    R1: {files['R1']['path']} ({files['R1']['size']} bytes)")
                if "R2" in files:
                    result.append(f"    R2: {files['R2']['path']} ({files['R2']['size']} bytes)")
        
        if single_files:
            result.append("\n单端测序文件：")
            for file_info in single_files:
                result.append(f"  {file_info['name']}: {file_info['path']} ({file_info['size']} bytes)")
        
        result.append(f"\n总计：{len(paired_files)} 个双端样本，{len(single_files)} 个单端文件")
        
        return "\n".join(result)
    
    except Exception as e:
        return f"查询FASTQ文件时发生错误：{str(e)}"

@tool(args_schema=GenomeQueryArgs)
def query_genome_info(genome_name: str, config_path: str = "config/genomes.json") -> str:
    """
    查询基因组配置信息
    
    遵循单一职责原则：专门处理基因组信息查询
    """
    try:
        if not os.path.exists(config_path):
            return f"错误：基因组配置文件 '{config_path}' 不存在"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        
        if genome_name not in genomes_config:
            available_genomes = list(genomes_config.keys())
            return f"错误：基因组 '{genome_name}' 不存在。可用基因组：{', '.join(available_genomes)}"
        
        genome_info = genomes_config[genome_name]
        
        result = [f"基因组 '{genome_name}' 信息："]
        result.append(f"  物种: {genome_info.get('species', '未知')}")
        result.append(f"  版本: {genome_info.get('version', '未知')}")
        result.append(f"  FASTA文件: {genome_info.get('fasta', '未配置')}")
        result.append(f"  GTF文件: {genome_info.get('gtf', '未配置')}")
        
        if 'fasta_url' in genome_info:
            result.append(f"  FASTA下载URL: {genome_info['fasta_url']}")
        if 'gtf_url' in genome_info:
            result.append(f"  GTF下载URL: {genome_info['gtf_url']}")
        
        # 检查本地文件是否存在
        fasta_path = genome_info.get('fasta', '')
        gtf_path = genome_info.get('gtf', '')
        
        if fasta_path and os.path.exists(fasta_path):
            result.append(f"  ✅ FASTA文件已存在")
        elif fasta_path:
            result.append(f"  ❌ FASTA文件不存在")
        
        if gtf_path and os.path.exists(gtf_path):
            result.append(f"  ✅ GTF文件已存在")
        elif gtf_path:
            result.append(f"  ❌ GTF文件不存在")
        
        return "\n".join(result)
    
    except json.JSONDecodeError:
        return f"错误：基因组配置文件 '{config_path}' 格式不正确"
    except Exception as e:
        return f"查询基因组信息时发生错误：{str(e)}"

# ============================================================================
# 配置管理工具组 - 遵循DRY原则
# ============================================================================

@tool(args_schema=NextflowConfigArgs)
def update_nextflow_param(param_name: str, param_value: Any) -> str:
    """
    更新单个nextflow参数
    
    应用KISS原则：简单的单参数更新
    """
    try:
        # 验证参数名称
        valid_params = [
            "srr_ids", "local_genome_path", "local_gtf_path", 
            "download_genome_url", "download_gtf_url", "local_fastq_files",
            "data", "run_download_srr", "run_download_genome", 
            "run_build_star_index", "run_fastp", "run_star_align", "run_featurecounts"
        ]
        
        if param_name not in valid_params:
            return f"错误：无效的参数名 '{param_name}'。有效参数：{', '.join(valid_params)}"
        
        # 类型验证
        if param_name.startswith("run_") and not isinstance(param_value, bool):
            return f"错误：参数 '{param_name}' 必须是布尔值"
        
        return f"✅ 参数 '{param_name}' 已更新为: {param_value}"
    
    except Exception as e:
        return f"更新参数时发生错误：{str(e)}"

@tool(args_schema=BatchConfigArgs)
def batch_update_nextflow_config(config_updates: Dict[str, Any]) -> str:
    """
    批量更新nextflow配置
    
    遵循DRY原则：复用单参数更新逻辑
    """
    try:
        results = []
        for param_name, param_value in config_updates.items():
            result = update_nextflow_param(param_name, param_value)
            results.append(result)
        
        success_count = sum(1 for r in results if r.startswith("✅"))
        error_count = len(results) - success_count
        
        summary = f"批量更新完成：{success_count} 个成功，{error_count} 个失败"
        return summary + "\n\n" + "\n".join(results)
    
    except Exception as e:
        return f"批量更新配置时发生错误：{str(e)}"

# ============================================================================
# 模式控制工具组 - 遵循开放封闭原则
# ============================================================================

@tool(args_schema=ModeSwitch)
def switch_to_plan_mode(target_mode: str, reason: str) -> str:
    """
    切换到计划模式
    
    遵循单一职责原则：专门处理模式切换
    """
    if target_mode != "plan":
        return f"错误：此工具只能切换到plan模式，收到：{target_mode}"
    
    return f"🔄 正在切换到计划模式...\n原因：{reason}\n✅ 模式切换成功！现在可以开始制定RNA-seq分析计划。"

@tool(args_schema=ModeSwitch)
def switch_to_execute_mode(target_mode: str, reason: str) -> str:
    """
    切换到执行模式
    
    遵循单一职责原则：专门处理模式切换
    """
    if target_mode != "execute":
        return f"错误：此工具只能切换到execute模式，收到：{target_mode}"
    
    return f"🔄 正在切换到执行模式...\n原因：{reason}\n✅ 模式切换成功！准备执行nextflow流程。"

# ============================================================================
# 执行控制工具组 - 遵循单一职责原则
# ============================================================================

@tool(args_schema=ExecutionArgs)
def execute_nextflow_pipeline(config_path: str = "config/nextflow.config", work_dir: str = "./work") -> str:
    """
    执行nextflow流程
    
    应用KISS原则：简单的流程执行
    """
    try:
        # 检查nextflow是否可用
        result = subprocess.run(["nextflow", "-version"], capture_output=True, text=True)
        if result.returncode != 0:
            return "错误：nextflow未安装或不可用"
        
        # 构建执行命令
        cmd = [
            "nextflow", "run", "main.nf",
            "-c", config_path,
            "-work-dir", work_dir
        ]
        
        return f"🚀 开始执行nextflow流程...\n命令：{' '.join(cmd)}\n⏳ 流程正在后台运行，请使用check_execution_status查看状态。"
    
    except Exception as e:
        return f"执行nextflow流程时发生错误：{str(e)}"

@tool
def check_execution_status() -> str:
    """
    检查执行状态
    
    应用KISS原则：简单的状态检查
    """
    try:
        # 检查是否有nextflow进程在运行
        result = subprocess.run(["pgrep", "-f", "nextflow"], capture_output=True, text=True)
        
        if result.returncode == 0:
            return "✅ Nextflow流程正在运行中..."
        else:
            return "⏹️ 没有检测到正在运行的nextflow流程"
    
    except Exception as e:
        return f"检查执行状态时发生错误：{str(e)}"

@tool
def get_current_nextflow_config() -> str:
    """
    获取当前nextflow配置状态
    
    遵循单一职责原则：专门获取配置信息
    """
    try:
        # 这里应该从状态管理中获取配置，暂时返回默认配置
        default_config = {
            "srr_ids": "",
            "local_genome_path": "",
            "local_gtf_path": "",
            "download_genome_url": "",
            "download_gtf_url": "",
            "local_fastq_files": "",
            "data": "./data",
            "run_download_srr": False,
            "run_download_genome": False,
            "run_build_star_index": False,
            "run_fastp": False,
            "run_star_align": False,
            "run_featurecounts": False
        }
        
        result = ["当前nextflow配置："]
        for key, value in default_config.items():
            result.append(f"  {key}: {value}")
        
        return "\n".join(result)
    
    except Exception as e:
        return f"获取配置时发生错误：{str(e)}"