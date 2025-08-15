import operator
from typing import Annotated, List, Dict, Any, Optional
from langchain_core.messages import AnyMessage
from typing_extensions import TypedDict

class AgentState(TypedDict):
    """
    RNA-seq分析Agent的核心状态管理
    
    遵循单一职责原则：每个字段承担明确的状态管理职责
    重构为支持Plan-Execute任务队列系统
    """
    # === 核心状态管理 ===
    mode: str  # "normal", "plan", "execute"
    messages: Annotated[List[AnyMessage], operator.add]
    
    # === Plan-Execute任务队列系统 ===
    task_queue: List[Dict[str, Any]]     # 待执行任务队列
    completed_tasks: List[Dict[str, Any]] # 已完成任务及结果
    current_task_index: int              # 当前执行任务索引
    
    # === 执行上下文管理 ===
    execution_context: Dict[str, Any]    # 执行上下文信息
    # 包含字段：
    # - iteration_count: 规划轮数
    # - tools_executed: 已执行工具列表  
    # - info_collected: 已收集信息状态
    # - config_completeness: 配置完整度百分比
    # - last_decisions: 最近的决策历史
    # - completion_criteria: 完成标准
    
    # === 原有字段保持兼容 ===
    current_stage: str               # 保持向后兼容
    tool_execution_history: List[str] # 保持向后兼容
    
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
    
    重构为支持Plan-Execute任务队列系统
    """
    return AgentState(
        mode="normal",
        messages=[],
        
        # === Plan-Execute任务队列系统初始化 ===
        task_queue=[],                   # 空任务队列
        completed_tasks=[],              # 空完成任务列表
        current_task_index=0,            # 从第一个任务开始
        
        # === 执行上下文初始化 ===
        execution_context={
            "iteration_count": 0,        # 规划轮数
            "tools_executed": [],        # 已执行工具
            "info_collected": {          # 信息收集状态
                "fastq_files": False,
                "genome_info": False,
                "basic_config": False
            },
            "config_completeness": 0,    # 配置完整度百分比
            "last_decisions": [],        # 最近决策历史
            "completion_criteria": {     # 完成标准
                "max_iterations": 5,
                "min_completeness": 80,
                "required_tools": ["query_fastq_files", "query_genome_info"],
                "required_info": ["fastq_files", "genome_info"]
            }
        },
        
        # === 保持向后兼容 ===
        current_stage="user_input",
        tool_execution_history=[],
        
        nextflow_config={
            # 默认nextflow参数
            "srr_ids": "",
            "local_genome_path": "",
            "local_gtf_path": "",
            "download_genome_url": "",
            "download_gtf_url": "",
            "local_fastq_files": "",
            "genome_version": "",
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

def update_current_stage(state: AgentState, new_stage: str, reason: str = "") -> AgentState:
    """
    更新当前执行阶段
    
    这是防止工具循环调用的核心机制
    """
    valid_stages = [
        "user_input", "planning_init", "tools_executing", 
        "tools_completed", "plan_ready", "executing_pipeline", 
        "execution_monitoring"
    ]
    
    if new_stage not in valid_stages:
        raise ValueError(f"Invalid stage: {new_stage}. Must be one of {valid_stages}")
    
    updated_state = state.copy()
    old_stage = updated_state.get("current_stage", "unknown")
    updated_state["current_stage"] = new_stage
    
    # 记录阶段转换日志
    if reason:
        log_message = f"阶段转换: {old_stage} → {new_stage} ({reason})"
    else:
        log_message = f"阶段转换: {old_stage} → {new_stage}"
    
    updated_state["execution_log"].append(f"[STAGE] {log_message}")
    
    return updated_state

def add_tool_to_execution_history(state: AgentState, tool_name: str) -> AgentState:
    """
    添加工具到执行历史，防止重复调用
    """
    updated_state = state.copy()
    if tool_name not in updated_state["tool_execution_history"]:
        updated_state["tool_execution_history"].append(tool_name)
    
    return updated_state

def has_tool_been_executed(state: AgentState, tool_name: str) -> bool:
    """
    检查工具是否已经被执行过
    """
    return tool_name in state.get("tool_execution_history", [])

def clear_tool_execution_history(state: AgentState) -> AgentState:
    """
    清空工具执行历史（用于新的计划周期）
    """
    updated_state = state.copy()
    updated_state["tool_execution_history"] = []
    return updated_state

def should_execute_tools(state: AgentState) -> bool:
    """
    判断当前阶段是否应该执行工具
    
    所有模式（normal/plan/execute）都支持工具调用
    """
    # 所有模式都允许工具调用
    return True

# ============================================================================
# Plan-Execute任务队列系统核心函数
# ============================================================================

def add_task_to_queue(state: AgentState, task: Dict[str, Any]) -> AgentState:
    """
    添加任务到队列
    """
    updated_state = state.copy()
    updated_state["task_queue"].append(task)
    return updated_state

def get_next_task(state: AgentState) -> Optional[Dict[str, Any]]:
    """
    获取下一个待执行任务
    """
    task_queue = state.get("task_queue", [])
    current_index = state.get("current_task_index", 0)
    
    if current_index < len(task_queue):
        return task_queue[current_index]
    return None

def mark_task_completed(state: AgentState, task_result: Dict[str, Any]) -> AgentState:
    """
    标记任务完成并记录结果
    """
    updated_state = state.copy()
    
    # 添加到已完成任务列表
    updated_state["completed_tasks"].append(task_result)
    
    # 更新任务索引
    updated_state["current_task_index"] += 1
    
    # 更新执行上下文
    context = updated_state["execution_context"]
    tool_name = task_result.get("tool_name", "unknown")
    if tool_name not in context["tools_executed"]:
        context["tools_executed"].append(tool_name)
    
    return updated_state

def is_task_queue_empty(state: AgentState) -> bool:
    """
    检查任务队列是否为空（所有任务已执行完成）
    """
    task_queue = state.get("task_queue", [])
    current_index = state.get("current_task_index", 0)
    
    return current_index >= len(task_queue)

def clear_task_queue(state: AgentState) -> AgentState:
    """
    清空任务队列（开始新的规划周期）
    """
    updated_state = state.copy()
    updated_state["task_queue"] = []
    updated_state["completed_tasks"] = []
    updated_state["current_task_index"] = 0
    return updated_state


def increment_iteration_count(state: AgentState) -> AgentState:
    """
    增加规划轮数计数
    """
    updated_state = state.copy()
    context = updated_state["execution_context"]
    context["iteration_count"] += 1
    
    # 记录本轮开始
    updated_state["execution_log"].append(f"[ITERATION] 开始第 {context['iteration_count']} 轮规划")
    
    return updated_state

def add_decision_to_history(state: AgentState, decision: str) -> AgentState:
    """
    添加决策到历史记录
    """
    updated_state = state.copy()
    context = updated_state["execution_context"]
    
    # 保持最近3个决策
    context["last_decisions"].append(decision)
    if len(context["last_decisions"]) > 3:
        context["last_decisions"] = context["last_decisions"][-3:]
    
    return updated_state

def calculate_config_completeness(state: AgentState) -> int:
    """
    计算配置完整度百分比
    """
    nextflow_config = state.get("nextflow_config", {})
    fastq_info = state.get("fastq_info", {})
    genome_info = state.get("genome_info", {})
    
    total_score = 0
    max_score = 100
    
    # FASTQ文件配置 (30分)
    if nextflow_config.get("local_fastq_files") or nextflow_config.get("srr_ids"):
        total_score += 30
    
    # 基因组配置 (30分)
    if (nextflow_config.get("local_genome_path") and nextflow_config.get("local_gtf_path")) or \
       nextflow_config.get("genome_version"):
        total_score += 30
    
    # 基本流程配置 (20分)
    process_params = ["run_fastp", "run_star_align", "run_featurecounts"]
    enabled_processes = sum(1 for param in process_params if nextflow_config.get(param))
    total_score += int(20 * enabled_processes / len(process_params))
    
    # 信息收集完成度 (20分)
    if fastq_info:
        total_score += 10
    if genome_info:
        total_score += 10
    
    return min(total_score, max_score)

def should_continue_planning(state: AgentState) -> tuple[bool, str]:
    """
    基于客观标准判断是否应该继续规划
    
    这是防死循环的核心函数
    """
    context = state.get("execution_context", {})
    criteria = context.get("completion_criteria", {})
    
    # 1. 检查最大迭代次数
    max_iterations = criteria.get("max_iterations", 5)
    current_iteration = context.get("iteration_count", 0)
    if current_iteration >= max_iterations:
        return False, f"达到最大迭代次数 ({max_iterations})"
    
    # 2. 检查配置完整度
    completeness = calculate_config_completeness(state)
    min_completeness = criteria.get("min_completeness", 80)
    if completeness >= min_completeness:
        return False, f"配置完整度足够 ({completeness}% >= {min_completeness}%)"
    
    # 3. 检查重复决策
    last_decisions = context.get("last_decisions", [])
    if len(last_decisions) >= 3 and len(set(last_decisions[-3:])) == 1:
        return False, f"检测到重复决策: {last_decisions[-1]}"
    
    # 4. 检查必要工具执行完成
    required_tools = criteria.get("required_tools", [])
    tools_executed = context.get("tools_executed", [])
    if all(tool in tools_executed for tool in required_tools):
        return False, "必要工具已全部执行完成"
    
    return True, "继续规划"

def update_config_completeness(state: AgentState) -> AgentState:
    """
    更新配置完整度
    """
    updated_state = state.copy()
    completeness = calculate_config_completeness(state)
    updated_state["execution_context"]["config_completeness"] = completeness
    return updated_state