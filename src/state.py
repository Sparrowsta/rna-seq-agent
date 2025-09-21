
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field
from .config.default_tool_params import DEFAULT_FASTP_PARAMS, DEFAULT_STAR_PARAMS, DEFAULT_FEATURECOUNTS_PARAMS, DEFAULT_HISAT2_PARAMS

# ==================== 统一Agent状态 ====================

class AgentState(BaseModel):
    """RNA-seq Agent统一状态 - 包含所有节点的字段"""
    
    # === 基础通信字段 ===
    messages: List[Any] = Field(default_factory=list, description="完整的LangGraph消息历史")
    input: str = Field(default="", description="用户当前输入")
    response: str = Field(default="", description="当前响应")
    status: str = Field(default="normal", description="系统状态: normal/plan/confirm/execute")
    
    # === Normal模式字段 ===
    routing_decision: str = Field(default="", description="路由决策: normal/execute/end (execute替代旧的plan)")
    query_response: str = Field(default="", description="查询响应内容")
    user_requirements: Dict[str, Any] = Field(default_factory=dict, description="从用户输入中提取的结构化配置需求")
    
    # === Detect检测字段 ===
    query_results: Dict[str, Any] = Field(default_factory=dict, description="系统检测结果")
    query_summary: str = Field(default="", description="检测结果总结")
    results_timestamp: str = Field(default="", description="结果目录时间戳")
    results_dir: str = Field(default="", description="完整结果目录路径")
    base_data_path: str = Field(default="", description="基础数据路径")
    
    # === Prepare配置字段 ===
    nextflow_config: Dict[str, Any] = Field(default_factory=dict, description="Nextflow配置参数")
    resource_config: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="各进程的CPU和内存资源配置")
    config_reasoning: str = Field(default="", description="配置决策理由")
    
    # === Modify修改字段 ===
    modification_history: List[Dict[str, Any]] = Field(default_factory=list, description="配置修改历史记录")
    # Modify入口临时携带的修改需求（来自Confirm节点）
    modify_requirements: Dict[str, Any] = Field(default_factory=dict, description="用户在/modify中输入的修改需求")
    
    # === UserConfirm确认字段 ===
    user_decision: str = Field(default="", description="用户决策: execute/replan/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")

    # === 执行模式（用于路由与展示） ===
    execution_mode: str = Field(default="single", description="执行模式: single/optimized/batch_optimize/yolo")
    current_step: str = Field(default="", description="当前执行步骤: fastp/star/featurecounts/analysis")
    completed_steps: List[str] = Field(default_factory=list, description="已完成的步骤列表")
    pipeline_progress: Dict[str, Any] = Field(default_factory=dict, description="流水线执行进度跟踪")


    # === 工具执行结果（保留这些字段，它们不是参数管理字段） ===
    fastp_results: Dict[str, Any] = Field(default_factory=dict, description="FastP执行产生的输出与摘要；遵循Prompt约定")
    star_results: Dict[str, Any] = Field(default_factory=dict, description="STAR执行产生的输出与摘要；遵循Prompt约定")  
    hisat2_results: Dict[str, Any] = Field(default_factory=dict, description="HISAT2执行产生的输出与摘要；遵循Prompt约定")
    featurecounts_results: Dict[str, Any] = Field(default_factory=dict, description="FeatureCounts执行产生的输出与摘要；遵循Prompt约定")

    # === 当前工具参数（仅保留当前参数，删除历史和版本管理） ===
    fastp_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FASTP_PARAMS.copy(), description="当前FastP运行参数")
    star_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_STAR_PARAMS.copy(), description="当前STAR比对参数")
    hisat2_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_HISAT2_PARAMS.copy(), description="当前HISAT2比对参数")
    featurecounts_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FEATURECOUNTS_PARAMS.copy(), description="当前FeatureCounts定量参数")

    # === 工具优化字段 ===
    fastp_optimization: Dict[str, Any] = Field(default_factory=dict, description="FastP优化结果")
    star_optimization: Dict[str, Any] = Field(default_factory=dict, description="STAR优化结果")
    hisat2_optimization: Dict[str, Any] = Field(default_factory=dict, description="HISAT2优化结果")
    featurecounts_optimization: Dict[str, Any] = Field(default_factory=dict, description="FeatureCounts优化结果")

    # === 工具优化参数变更字段（供跨节点引用优化历史） ===
    fastp_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="FastP优化参数变更记录（仅包含改动项）")
    star_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="STAR优化参数变更记录（仅包含改动项）")
    hisat2_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="HISAT2优化参数变更记录（仅包含改动项）")
    featurecounts_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="FeatureCounts优化参数变更记录（仅包含改动项）")

    # === 工具优化建议字段 ===
    fastp_optimization_suggestions: str = Field(default="", description="FastP优化建议")
    star_optimization_suggestions: str = Field(default="", description="STAR优化建议")
    hisat2_optimization_suggestions: str = Field(default="", description="HISAT2优化建议")
    featurecounts_optimization_suggestions: str = Field(default="", description="FeatureCounts优化建议")
    
    # === 工具优化历史字段（支持多次历史保存） ===
    fastp_optimization_history: List[Dict[str, Any]] = Field(default_factory=list, description="FastP优化参数历史列表，保存最近5次优化记录")
    star_optimization_history: List[Dict[str, Any]] = Field(default_factory=list, description="STAR优化参数历史列表，保存最近5次优化记录")
    hisat2_optimization_history: List[Dict[str, Any]] = Field(default_factory=list, description="HISAT2优化参数历史列表，保存最近5次优化记录")
    featurecounts_optimization_history: List[Dict[str, Any]] = Field(default_factory=list, description="FeatureCounts优化参数历史列表，保存最近5次优化记录")
    # === 工作流集成字段 ===
    workflow_status: str = Field(default="", description="整体流程状态")
    workflow_checkpoints: Dict[str, Any] = Field(default_factory=dict, description="检查点数据")
    optimization_suggestions: Dict[str, Any] = Field(default_factory=dict, description="优化建议")
    user_param_overrides: Dict[str, Any] = Field(default_factory=dict, description="用户参数覆盖")
    workflow_summary: str = Field(default="", description="流程总结报告")
    
    # === Analysis分析字段 ===
    overall_summary: str = Field(default="", description="流水线执行和数据质量的整体摘要，包括成功状态和完成情况")
    key_findings: List[str] = Field(default_factory=list, description="基于数据分析的关键发现和模式，重要的数据洞察和生物学意义")
    sample_health_assessment: str = Field(default="", description="各样本的健康度评估和问题标记，包括PASS/WARN/FAIL状态判断")
    quality_metrics_analysis: str = Field(default="", description="FastP、比对、定量等步骤的质量指标专业解读和数据模式分析")
    optimization_recommendations: List[str] = Field(default_factory=list, description="具体的参数调优和实验改进建议，基于数据质量的可行建议")
    risk_warnings: List[str] = Field(default_factory=list, description="数据使用和后续分析的注意事项，潜在风险和限制因素")
    next_steps: List[str] = Field(default_factory=list, description="建议的后续分析方向和步骤，包括差异分析、功能富集等")


# ==================== 节点响应模型 - 用于特定节点的结构化输出 ====================


class NormalResponse(BaseModel):
    """Normal节点的精简响应格式 - 兼容create_react_agent工具响应"""
    routing_decision: str = Field(default="", description="路由决策: normal/plan/end")
    query_response: str = Field(description="工具调用的完整结果")
    user_requirements: Dict[str, Any] = Field(default_factory=dict, description="从用户输入中提取的结构化配置需求")

class PrepareResponse(BaseModel):
    """Prepare节点的精简响应格式"""
    nextflow_config: Dict[str, Any] = Field(default_factory=dict, description="生成的Nextflow配置参数")
    resource_config: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="各进程的CPU和内存资源配置")
    config_reasoning: str = Field(default="", description="配置决策理由")

class FastpResponse(BaseModel):
    """FastP节点的Agent响应格式"""
    fastp_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FASTP_PARAMS.copy(), description="优化后的FastP参数")
    fastp_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    fastp_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="仅包含改变了的优化参数")
    fastp_results: Dict[str, Any] = Field(default_factory=dict, description="FastP执行产生的输出与摘要；遵循Prompt约定")

class StarResponse(BaseModel):
    """STAR节点的Agent响应格式"""
    star_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_STAR_PARAMS.copy(), description="优化后的STAR参数")
    star_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    star_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="仅包含改变了的优化参数")
    star_results: Dict[str, Any] = Field(default_factory=dict, description="STAR执行产生的输出与摘要；遵循Prompt约定")

class FeaturecountsResponse(BaseModel):
    """FeatureCounts节点的Agent响应格式"""
    featurecounts_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_FEATURECOUNTS_PARAMS.copy(), description="优化后的FeatureCounts参数")
    featurecounts_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    featurecounts_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="仅包含改变了的优化参数")
    featurecounts_results: Dict[str, Any] = Field(default_factory=dict, description="FeatureCounts执行产生的输出与摘要；遵循Prompt约定")

class Hisat2Response(BaseModel):
    """HISAT2节点的Agent响应格式"""
    hisat2_params: Dict[str, Any] = Field(default_factory=lambda: DEFAULT_HISAT2_PARAMS.copy(), description="优化后的HISAT2参数")
    hisat2_optimization_suggestions: str = Field(description="详细的优化理由、数据支撑和预期效果")
    hisat2_optimization_params: Dict[str, Any] = Field(default_factory=dict, description="仅包含改变了的优化参数")
    hisat2_results: Dict[str, Any] = Field(default_factory=dict, description="HISAT2执行产生的输出与摘要；遵循Prompt约定")

class AnalysisResponse(BaseModel):
    """RNA-seq分析结果的结构化响应格式"""
    overall_summary: str = Field(default="", description="流水线执行和数据质量的整体摘要，包括成功状态和完成情况")
    key_findings: List[str] = Field(default_factory=list, description="基于数据分析的关键发现和模式，重要的数据洞察和生物学意义")
    sample_health_assessment: str = Field(default="", description="各样本的健康度评估和问题标记，包括PASS/WARN/FAIL状态判断")
    quality_metrics_analysis: str = Field(default="", description="FastP、比对、定量等步骤的质量指标专业解读和数据模式分析")
    optimization_recommendations: List[str] = Field(default_factory=list, description="具体的参数调优和实验改进建议，基于数据质量的可行建议")
    risk_warnings: List[str] = Field(default_factory=list, description="数据使用和后续分析的注意事项，潜在风险和限制因素")
    next_steps: List[str] = Field(default_factory=list, description="建议的后续分析方向和步骤，包括差异分析、功能富集等")
