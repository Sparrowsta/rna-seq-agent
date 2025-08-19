
from typing import List, Dict, Any
from pydantic import BaseModel, Field

# ==================== 基础Pydantic状态类 ====================

class BaseAgentState(BaseModel):
    """基础Agent状态 - 包含所有节点共享的公共字段"""
    
    # 公共字段 - 所有节点都需要
    messages: List[Dict[str, str]] = Field(default=[], description="对话历史，格式: [{'role': 'user', 'content': '...'}]")
    response: str = Field(default="", description="当前阶段响应")
    status: str = Field(default="pending", description="当前状态: pending, planning, detecting, preparing, ready, replanning, confirming, executing")

class NormalNodeState(BaseAgentState):
    """Normal节点状态 - 用户交互和信息查询"""
    
    # Normal节点专用字段
    query_type: str = Field(default="", description="查询类型: info/analysis/help/plan")
    routing_decision: str = Field(default="", description="路由决策: normal/plan")
    query_response: str = Field(default="", description="查询响应内容")
    user_intent: str = Field(default="", description="用户意图解析")
    suggested_actions: List[str] = Field(default=[], description="建议的操作列表")

class UserCommunicationNodeState(BaseAgentState):
    """User Communication节点状态 - 用户交互和功能执行"""
    
    # 继承Normal节点输出
    query_type: str = Field(default="", description="查询类型: info/analysis/help/plan")
    routing_decision: str = Field(default="", description="路由决策: normal/plan")
    query_response: str = Field(default="", description="查询响应内容")
    user_intent: str = Field(default="", description="用户意图解析")
    suggested_actions: List[str] = Field(default=[], description="建议的操作列表")
    
    # User Communication节点专用字段
    execution_result: Dict[str, Any] = Field(default={}, description="功能执行结果")
    user_feedback: str = Field(default="", description="用户反馈")
    next_action: str = Field(default="", description="下一步操作建议")

class PlanNodeState(BaseAgentState):
    """Plan节点状态 - 需求理解和计划制定"""
    
    # Plan节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")

class DetectNodeState(BaseAgentState):
    """Detect节点状态 - 系统信息检测"""
    
    # 继承Plan节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    
    # Detect节点输出
    query_results: Dict[str, Any] = Field(default={}, description="检测到的系统信息")
    query_summary: str = Field(default="", description="检测结果总结")

class PrepareNodeState(BaseAgentState):
    """Prepare节点状态 - 配置准备和参数生成"""
    
    # 继承前面节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    query_results: Dict[str, Any] = Field(default={}, description="检测到的系统信息")
    query_summary: str = Field(default="", description="检测结果总结")
    
    # Prepare节点输出
    nextflow_config: Dict[str, Any] = Field(default={}, description="最终Nextflow配置")
    config_reasoning: str = Field(default="", description="配置决策理由")

class ReplanNodeState(BaseAgentState):
    """Replan节点状态 - 用户修改需求分析和路由决策"""
    
    # 继承前面节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    query_results: Dict[str, Any] = Field(default={}, description="检测到的系统信息")
    query_summary: str = Field(default="", description="检测结果总结")
    nextflow_config: Dict[str, Any] = Field(default={}, description="Nextflow配置")
    
    # Replan节点输出
    user_modification_input: str = Field(default="", description="用户修改请求")
    modification_intent: Dict[str, Any] = Field(default={}, description="解析的修改意图")
    modification_mode: str = Field(default="", description="修改模式: incremental/full")
    routing_decision: str = Field(default="", description="路由决策: detect")
    routing_reason: str = Field(default="", description="路由理由")

class UserConfirmState(BaseAgentState):
    """用户确认节点状态 - 展示配置并等待用户决策"""
    
    # 继承前面所有节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    query_results: Dict[str, Any] = Field(default={}, description="检测到的系统信息")
    query_summary: str = Field(default="", description="检测结果总结")
    nextflow_config: Dict[str, Any] = Field(default={}, description="Nextflow配置")
    config_reasoning: str = Field(default="", description="配置决策理由")
    
    # UserConfirm节点输出
    user_decision: str = Field(default="", description="用户决策: execute/modify/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")

class ExecuteNodeState(BaseAgentState):
    """执行节点状态 - Nextflow命令构建和执行"""
    
    # 继承前面所有节点输出
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    query_results: Dict[str, Any] = Field(default={}, description="检测到的系统信息")
    query_summary: str = Field(default="", description="检测结果总结")
    nextflow_config: Dict[str, Any] = Field(default={}, description="Nextflow配置")
    config_reasoning: str = Field(default="", description="配置决策理由")
    user_decision: str = Field(default="", description="用户决策")
    
    # Execute节点输出
    nextflow_command: str = Field(default="", description="构建的nextflow命令")
    execution_status: str = Field(default="", description="执行状态: building/running/completed/failed")
    execution_output: str = Field(default="", description="执行输出日志")
    execution_result: Dict[str, Any] = Field(default={}, description="执行结果摘要")
