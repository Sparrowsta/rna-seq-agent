
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field
from .config.default_tool_params import DEFAULT_FASTP_PARAMS, DEFAULT_STAR_PARAMS, DEFAULT_FEATURECOUNTS_PARAMS

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
    # Modify入口临时携带的修改需求（来自Confirm节点）
    modify_requirements: Dict[str, Any] = Field(default={}, description="用户在/modify中输入的修改需求")
    
    # === UserConfirm确认字段 ===
    user_decision: str = Field(default="", description="用户决策: execute/replan/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")

    # === 执行模式（用于路由与展示） ===
    execution_mode: str = Field(default="single", description="执行模式: single/optimized/batch_optimize/yolo")
    current_step: str = Field(default="", description="当前执行步骤: fastp/star/featurecounts/analysis")
    completed_steps: List[str] = Field(default=[], description="已完成的步骤列表")
    pipeline_progress: Dict[str, Any] = Field(default={}, description="流水线执行进度跟踪")

    # === 批次优化模式字段 ===
    batch_optimization_mode: bool = Field(default=False, description="是否启用批次优化模式 - 当execution_mode为batch_optimize时自动设为True")
    batch_optimizations: Dict[str, Dict[str, Any]] = Field(default={}, description="批次优化收集的所有工具优化参数")
    batch_optimization_round: int = Field(default=1, description="批次优化轮次")
    batch_optimization_complete: bool = Field(default=False, description="批次优化是否完成")

    # === 工具执行结果（保留这些字段，它们不是参数管理字段） ===
    fastp_results: Dict[str, Any] = Field(default={}, description="FastP执行结果数据")
    star_results: Dict[str, Any] = Field(default={}, description="STAR比对结果数据")  
    featurecounts_results: Dict[str, Any] = Field(default={}, description="FeatureCounts定量结果数据")

    # === 当前工具参数（仅保留当前参数，删除历史和版本管理） ===
    fastp_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FASTP_PARAMS.copy(), description="当前FastP运行参数")
    star_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_STAR_PARAMS.copy(), description="当前STAR比对参数")
    featurecounts_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FEATURECOUNTS_PARAMS.copy(), description="当前FeatureCounts定量参数")

    # === 工具优化字段 ===
    fastp_optimization: Dict[str, Any] = Field(default={}, description="FastP优化结果")
    star_optimization: Dict[str, Any] = Field(default={}, description="STAR优化结果")
    featurecounts_optimization: Dict[str, Any] = Field(default={}, description="FeatureCounts优化结果")

    # === 工具优化参数变更字段（供跨节点引用优化历史） ===
    fastp_optimization_params: Dict[str, Any] = Field(default={}, description="FastP优化参数变更记录（仅包含改动项）")
    star_optimization_params: Dict[str, Any] = Field(default={}, description="STAR优化参数变更记录（仅包含改动项）")
    featurecounts_optimization_params: Dict[str, Any] = Field(default={}, description="FeatureCounts优化参数变更记录（仅包含改动项）")

    # === 工具优化建议字段 ===
    fastp_optimization_suggestions: str = Field(default="", description="FastP优化建议")
    star_optimization_suggestions: str = Field(default="", description="STAR优化建议")
    featurecounts_optimization_suggestions: str = Field(default="", description="FeatureCounts优化建议")
    
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

# LLM 智能分析结构化输出模型
class LLMAnalysisModel(BaseModel):
    """LLM 智能分析的结构化输出模型"""
    global_summary: str  # 3-5句面向非技术读者的总结
    key_findings: List[str]  # 每条包含具体数据的关键发现
    per_sample_flags: List[Dict[str, Any]]  # sample_id, issues, severity
    recommendations: List[Dict[str, str]]  # type, title, detail
    risks: List[str]  # 潜在风险与注意事项
    report_md: Optional[str] = None  # 可选的额外Markdown片段
    debug_notes: Optional[List[str]] = None  # 调试模式时的补充信息

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
    """FastP节点的Agent响应格式"""
    fastp_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FASTP_PARAMS.copy(), description="优化后的FastP参数")
    fastp_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    fastp_optimization_params: Dict[str, Any] = Field(default={}, description="仅包含改变了的优化参数")
    results: Dict[str, Any] = Field(default={}, description="FastP执行产生的输出与摘要；遵循Prompt约定")

class StarResponse(BaseModel):
    """STAR节点的Agent响应格式"""
    star_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_STAR_PARAMS.copy(), description="优化后的STAR参数")
    star_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    star_optimization_params: Dict[str, Any] = Field(default={}, description="仅包含改变了的优化参数")
    results: Dict[str, Any] = Field(default={}, description="STAR执行产生的输出与摘要；遵循Prompt约定")

class FeaturecountsResponse(BaseModel):
    """FeatureCounts节点的Agent响应格式"""
    featurecounts_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FEATURECOUNTS_PARAMS.copy(), description="优化后的FeatureCounts参数")
    featurecounts_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    featurecounts_optimization_params: Dict[str, Any] = Field(default={}, description="仅包含改变了的优化参数")
    results: Dict[str, Any] = Field(default={}, description="FeatureCounts执行产生的输出与摘要；遵循Prompt约定")
