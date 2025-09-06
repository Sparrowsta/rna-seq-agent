
from typing import List, Dict, Any, Union
from pydantic import BaseModel, Field

# ==================== 统一Agent状态 ====================

class AgentState(BaseModel):
    """RNA-seq Agent统一状态 - 包含所有节点的字段"""
    
    # === 基础通信字段 ===
    messages: List[Any] = Field(default=[], description="完整的LangChain消息历史")
    input: str = Field(default="", description="用户当前输入")
    response: str = Field(default="", description="当前响应")
    status: str = Field(default="normal", description="系统状态: normal/plan/confirm/execute")
    
    # === Normal模式字段 ===
    routing_decision: str = Field(default="", description="路由决策: normal/plan/end")
    query_response: str = Field(default="", description="查询响应内容")
    
    # === Plan分析字段 ===
    plan: List[List[str]] = Field(default=[], description="并行任务组列表，每个组内任务串行执行")
    group_descriptions: List[str] = Field(default=[], description="任务组功能描述")
    execution_strategy: str = Field(default="parallel", description="执行策略")
    analysis_intent: str = Field(default="", description="分析目标意图")
    user_requirements: Dict[str, Any] = Field(default={}, description="初始用户配置需求(从user_communication来)")
    modify_requirements: Dict[str, Any] = Field(default={}, description="用户修改配置需求(从user_confirm的/modify来)")
    
    # === Detect检测字段 ===
    query_results: Dict[str, Any] = Field(default={}, description="系统检测结果")
    query_summary: str = Field(default="", description="检测结果总结")
    results_timestamp: str = Field(default="", description="结果目录时间戳")
    results_dir: str = Field(default="", description="完整结果目录路径")
    base_data_path: str = Field(default="", description="基础数据路径")
    
    # === Prepare配置字段 ===
    nextflow_config: Dict[str, Any] = Field(default={}, description="Nextflow配置参数")
    resource_config: Dict[str, Dict[str, Any]] = Field(default={}, description="各进程的CPU和内存资源配置")
    config_reasoning: str = Field(default="", description="配置决策理由")
    
    # === UserConfirm确认字段 ===
    user_decision: str = Field(default="", description="用户决策: execute/replan/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")

    # === 执行模式（用于路由与展示） ===
    execution_mode: str = Field(default="single", description="执行模式: single/optimized")

    # === 工具参数（对比展示用） ===
    # 目前按需支持 fastp；后续可扩展 star/featurecounts 同名字段
    fastp_default_params: Dict[str, Any] = Field(default={}, description="fastp 默认参数集（展示用）")
    fastp_optimized_params: Dict[str, Any] = Field(default={}, description="fastp 优化后的参数集（展示用）")
    fastp_current_params: Dict[str, Any] = Field(default={}, description="fastp 当前运行参数（稳定基线）")
    
    # === 参数版本管理 ===
    fastp_version: int = Field(default=1, description="fastp 参数版本号")
    fastp_version_history: List[Dict[str, Any]] = Field(default=[], description="fastp 参数历史版本记录")
    # 执行期参数对比辅助：用于确认面板内联展示 old -> new
    fastp_prev_params: Dict[str, Any] = Field(default={}, description="fastp 本次执行前使用的参数")
    fastp_applied_updates: Dict[str, Any] = Field(default={}, description="fastp 本次执行应用的参数差异（旧->新）")

# ==================== 子状态模型 - 用于特定节点的结构化输出 ====================

class NormalResponse(BaseModel):
    """Normal节点的精简响应格式 - 兼容create_react_agent工具响应"""
    query_response: str = Field(description="工具调用的完整结果")
    user_requirements: Dict[str, Any] = Field(default={}, description="从用户输入中提取的结构化配置需求")

class DetectResponse(BaseModel):
    """Detect节点的精简响应格式"""
    query_results: Dict[str, Any] = Field(default={}, description="系统检测结果")
    query_summary: str = Field(default="", description="检测结果总结")

class PrepareResponse(BaseModel):
    """Prepare节点的精简响应格式"""
    nextflow_config: Dict[str, Any] = Field(default={}, description="生成的Nextflow配置参数")
    resource_config: Dict[str, Dict[str, Any]] = Field(default={}, description="各进程的CPU和内存资源配置")
    config_reasoning: str = Field(default="", description="配置决策理由")

# （分析节点已移除，不再定义 AnalysisResponse）
