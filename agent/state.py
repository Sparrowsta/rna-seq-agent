import operator
from typing import Annotated, List, Dict, Any, Optional
from langchain_core.messages import AnyMessage
from typing_extensions import TypedDict

class AgentState(TypedDict):
    """
    RNA-seq分析Agent的核心状态管理
    
    遵循单一职责原则：每个字段承担明确的状态管理职责
    """
    # === 核心状态管理 ===
    mode: str  # "normal", "plan", "execute"
    messages: Annotated[List[AnyMessage], operator.add]
    
    # === 配置管理 ===
    nextflow_config: Dict[str, Any]  # nextflow流程参数
    genome_config: Dict[str, Any]    # 基因组配置信息
    
    # === 计划管理 ===
    plan: List[str]                  # 分析计划步骤列表
    current_step: int                # 当前执行的步骤索引
    plan_status: str                 # "draft", "confirmed", "executing", "completed"
    
    # === 数据信息 ===
    fastq_info: Dict[str, Any]       # FASTQ文件信息和元数据
    genome_info: Dict[str, Any]      # 选定基因组的详细信息
    
    # === 执行状态 ===
    execution_results: Dict[str, Any] # 各步骤的执行结果
    execution_status: str            # "idle", "running", "completed", "failed"
    execution_log: List[str]         # 执行日志记录

def create_initial_state() -> AgentState:
    """
    创建初始状态
    
    应用KISS原则：提供简单的初始化方法
    """
    return AgentState(
        mode="normal",
        messages=[],
        nextflow_config={
            # 默认nextflow参数
            "srr_ids": "",
            "local_genome_path": "",
            "local_gtf_path": "",
            "download_genome_url": "",
            "download_gtf_url": "",
            "local_fastq_files": "",
            "genome_version": "",  # 添加genome_version字段，初始为空
            "data": "./data",
            # 流程控制参数
            "run_download_srr": False,
            "run_download_genome": False,
            "run_build_star_index": False,
            "run_fastp": False,
            "run_star_align": False,
            "run_featurecounts": False
        },
        genome_config={},
        plan=[],
        current_step=0,
        plan_status="draft",
        fastq_info={},
        genome_info={},
        execution_results={},
        execution_status="idle",
        execution_log=[]
    )

def update_state_mode(state: AgentState, new_mode: str) -> AgentState:
    """
    更新状态模式
    
    遵循开放封闭原则：通过函数扩展功能而不修改状态结构
    """
    valid_modes = ["normal", "plan", "execute"]
    if new_mode not in valid_modes:
        raise ValueError(f"Invalid mode: {new_mode}. Must be one of {valid_modes}")
    
    updated_state = state.copy()
    updated_state["mode"] = new_mode
    return updated_state

def update_nextflow_config(state: AgentState, config_updates: Dict[str, Any]) -> AgentState:
    """
    更新nextflow配置
    
    应用DRY原则：统一的配置更新逻辑
    """
    updated_state = state.copy()
    updated_state["nextflow_config"].update(config_updates)
    return updated_state

def add_plan_step(state: AgentState, step: str) -> AgentState:
    """
    添加计划步骤
    
    遵循单一职责原则：专门处理计划管理
    """
    updated_state = state.copy()
    updated_state["plan"].append(step)
    return updated_state

def update_execution_status(state: AgentState, status: str, log_message: str = "") -> AgentState:
    """
    更新执行状态
    
    应用KISS原则：简单直接的状态更新
    """
    valid_statuses = ["idle", "running", "completed", "failed"]
    if status not in valid_statuses:
        raise ValueError(f"Invalid status: {status}. Must be one of {valid_statuses}")
    
    updated_state = state.copy()
    updated_state["execution_status"] = status
    
    if log_message:
        updated_state["execution_log"].append(f"[{status.upper()}] {log_message}")
    
    return updated_state