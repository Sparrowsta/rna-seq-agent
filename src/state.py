
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
    plan: List[str] = Field(default=[], description="分析步骤计划")
    analysis_intent: str = Field(default="", description="分析目标意图")
    user_requirements: Dict[str, Any] = Field(default={}, description="初始用户配置需求(从user_communication来)")
    replan_requirements: Dict[str, Any] = Field(default={}, description="重新规划配置需求(从user_confirm来)")
    
    # === Detect检测字段 ===
    query_results: Dict[str, Any] = Field(default={}, description="系统检测结果")
    query_summary: str = Field(default="", description="检测结果总结")
    
    # === Prepare配置字段 ===
    nextflow_config: Dict[str, Any] = Field(default={}, description="Nextflow配置参数")
    config_reasoning: str = Field(default="", description="配置决策理由")
    
    # === UserConfirm确认字段 ===
    user_decision: str = Field(default="", description="用户决策: execute/replan/cancel")
    confirmation_message: str = Field(default="", description="确认界面展示信息")
    
    # === Execute执行字段 ===
    nextflow_command: str = Field(default="", description="构建的Nextflow命令")
    execution_status: str = Field(default="", description="执行状态: building/running/completed/failed")
    execution_output: str = Field(default="", description="执行输出日志")
    execution_result: Dict[str, Any] = Field(default={}, description="执行结果摘要")

# ==================== 子状态模型 - 用于特定节点的结构化输出 ====================

class NormalResponse(BaseModel):
    """Normal节点的精简响应格式 - 兼容create_react_agent工具响应"""
    query_response: str = Field(description="工具调用的完整结果")
    user_requirements: Dict[str, Any] = Field(default={}, description="从用户输入中提取的结构化配置需求")

class PlanResponse(BaseModel):
    """Plan节点的精简响应格式"""
    plan: List[str] = Field(default=[], description="分析步骤计划")

class DetectResponse(BaseModel):
    """Detect节点的精简响应格式"""
    query_results: Dict[str, Any] = Field(default={}, description="系统检测结果")
    query_summary: str = Field(default="", description="检测结果总结")

class PrepareResponse(BaseModel):
    """Prepare节点的精简响应格式"""
    nextflow_config: Dict[str, Any] = Field(default={}, description="生成的Nextflow配置参数")
    config_reasoning: str = Field(default="", description="配置决策理由")
