
from typing import List, Dict, Any, Union, Optional, Literal
from pydantic import BaseModel, Field
from datetime import datetime

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
    user_requirements: Dict[str, Any] = Field(default={}, description="从用户输入中提取的结构化配置需求")
    
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
    
    # === Modify修改字段 ===
    modification_history: List[Dict[str, Any]] = Field(default=[], description="配置修改历史记录")
    
    # === UserConfirm确认字段 ===
    user_decision: str = Field(default="", description="用户决策: execute/replan/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")

    # === 执行模式（用于路由与展示） ===
    execution_mode: str = Field(default="single", description="执行模式: single/optimized")
    current_step: str = Field(default="", description="当前执行步骤: fastp/star/featurecounts/analysis")
    completed_steps: List[str] = Field(default=[], description="已完成的步骤列表")
    pipeline_progress: Dict[str, Any] = Field(default={}, description="流水线执行进度跟踪")

    # === FastP参数管理 ===
    fastp_params: Dict[str, Any] = Field(default={}, description="FastP运行参数")
    fastp_params_history: List[Dict[str, Any]] = Field(default=[], description="FastP参数执行历史")
    fastp_results: Dict[str, Any] = Field(default={}, description="FastP执行结果数据")
    fastp_version: int = Field(default=1, description="FastP参数版本号")
    fastp_version_history: List[Dict[str, Any]] = Field(default=[], description="FastP版本历史")
    fastp_optimization_suggestions: str = Field(default="", description="FastP参数优化建议")
    fastp_optimization_params: Dict[str, Any] = Field(default={}, description="FastP优化参数字典")
    
    # === STAR参数管理 ===
    star_params: Dict[str, Any] = Field(default={}, description="STAR比对参数")
    star_params_history: List[Dict[str, Any]] = Field(default=[], description="STAR参数执行历史")
    star_results: Dict[str, Any] = Field(default={}, description="STAR比对结果数据")
    star_version: int = Field(default=1, description="STAR参数版本号")
    star_version_history: List[Dict[str, Any]] = Field(default=[], description="STAR版本历史")
    star_optimization_suggestions: str = Field(default="", description="STAR参数优化建议")
    star_optimization_params: Dict[str, Any] = Field(default={}, description="STAR优化参数字典")
    
    # === FeatureCounts参数管理 ===
    featurecounts_params: Dict[str, Any] = Field(default={}, description="FeatureCounts定量参数")
    featurecounts_params_history: List[Dict[str, Any]] = Field(default=[], description="FeatureCounts参数执行历史")
    featurecounts_results: Dict[str, Any] = Field(default={}, description="FeatureCounts定量结果数据")
    featurecounts_version: int = Field(default=1, description="FeatureCounts参数版本号")
    featurecounts_version_history: List[Dict[str, Any]] = Field(default=[], description="FeatureCounts版本历史")
    featurecounts_optimization_suggestions: str = Field(default="", description="FeatureCounts参数优化建议")
    featurecounts_optimization_params: Dict[str, Any] = Field(default={}, description="FeatureCounts优化参数字典")

    # === 工作流集成字段 ===
    workflow_status: str = Field(default="", description="整体流程状态")
    workflow_checkpoints: Dict[str, Any] = Field(default={}, description="检查点数据")
    optimization_suggestions: Dict[str, Any] = Field(default={}, description="优化建议")
    user_param_overrides: Dict[str, Any] = Field(default={}, description="用户参数覆盖")
    workflow_summary: str = Field(default="", description="流程总结报告")
    
    # === Analysis分析字段 ===
    analysis_report: Dict[str, Any] = Field(default={}, description="综合分析报告数据")
    analysis_report_path: str = Field(default="", description="分析报告文件路径")
    workflow_statistics: Dict[str, Any] = Field(default={}, description="工作流统计数据")
    rna_seq_complete: bool = Field(default=False, description="RNA-seq流程完整标志")

# ==================== 节点响应模型 - 用于特定节点的结构化输出 ====================

class NormalResponse(BaseModel):
    """Normal节点的精简响应格式 - 兼容create_react_agent工具响应"""
    routing_decision: str = Field(default="", description="路由决策: normal/plan/end")
    query_response: str = Field(description="工具调用的完整结果")
    user_requirements: Dict[str, Any] = Field(default={}, description="从用户输入中提取的结构化配置需求")

class PrepareResponse(BaseModel):
    """Prepare节点的精简响应格式"""
    nextflow_config: Dict[str, Any] = Field(default={}, description="生成的Nextflow配置参数")
    resource_config: Dict[str, Dict[str, Any]] = Field(default={}, description="各进程的CPU和内存资源配置")
    config_reasoning: str = Field(default="", description="配置决策理由")

class FastpResponse(BaseModel):
    """FastP节点的响应格式"""
    status: str = Field(description="执行状态")
    summary: str = Field(description="执行总结")
    response: str = Field(description="节点响应消息")

class StarResponse(BaseModel):
    """STAR节点的响应格式"""
    status: str = Field(description="执行状态")
    summary: str = Field(description="执行总结")
    response: str = Field(description="节点响应消息")

class FeaturecountsResponse(BaseModel):
    """FeatureCounts节点的响应格式"""
    status: str = Field(description="执行状态")
    summary: str = Field(description="执行总结")
    response: str = Field(description="节点响应消息")