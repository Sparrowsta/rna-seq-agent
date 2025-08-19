
# ==================== 基础状态类 ====================

class BaseAgentState(dict):
    """基础Agent状态 - 包含所有节点共享的公共字段"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 公共字段 - 所有节点都需要
            "input": "",          # 用户原始输入
            "messages": [],       # 对话历史
            "response": "",       # 当前阶段响应
            "status": "pending"   # 当前状态: pending, planning, detecting, preparing, ready, replanning, confirming, executing
        })

# ==================== 节点专用状态类 ====================

class NormalNodeState(BaseAgentState):
    """Normal节点状态 - 用户交互和信息查询"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # Normal节点输出
            "query_type": "",            # 查询类型: info/analysis/help
            "routing_decision": "",      # 路由决策: normal/plan
            "query_response": "",        # 查询响应内容
        })

class PlanNodeState(BaseAgentState):
    """Plan节点状态 - 需求理解和计划制定"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # Plan节点输出
            "plan": [],              # 分析步骤计划
            "analysis_intent": "",   # 分析目标意图
        })

class DetectNodeState(BaseAgentState):
    """Detect节点状态 - 系统信息检测"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 继承Plan节点输出
            "plan": [],
            "analysis_intent": "",
            
            # Detect节点输出
            "query_results": {},     # 检测到的系统信息
            "query_summary": "",     # 检测结果总结
        })

class PrepareNodeState(BaseAgentState):
    """Prepare节点状态 - 配置准备和参数生成"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 继承前面节点输出
            "plan": [],
            "analysis_intent": "",
            "query_results": {},
            "query_summary": "",
            
            # Prepare节点输出
            "nextflow_config": {},   # 最终Nextflow配置
            "config_reasoning": "",  # 配置决策理由
        })

class ReplanNodeState(BaseAgentState):
    """Replan节点状态 - 用户修改需求分析和路由决策"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 继承前面节点输出
            "plan": [],
            "analysis_intent": "",
            "query_results": {},
            "query_summary": "",
            "nextflow_config": {},
            
            # Replan节点输出
            "user_modification_input": "",  # 用户修改请求
            "modification_intent": {},      # 解析的修改意图
            "modification_mode": "",        # 修改模式: incremental/full
            "routing_decision": "",         # 路由决策: detect
            "routing_reason": "",           # 路由理由
        })

class UserConfirmState(BaseAgentState):
    """用户确认节点状态 - 展示配置并等待用户决策"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 继承前面所有节点输出
            "plan": [],
            "analysis_intent": "",
            "query_results": {},
            "query_summary": "",
            "nextflow_config": {},
            "config_reasoning": "",
            
            # UserConfirm节点输出
            "user_decision": "",        # 用户决策: execute/modify/cancel
            "confirmation_message": "", # 确认界面展示信息
        })

class ExecuteNodeState(BaseAgentState):
    """执行节点状态 - Nextflow命令构建和执行"""
    
    def __init__(self):
        super().__init__()
        self.update({
            # 继承前面所有节点输出
            "plan": [],
            "analysis_intent": "",
            "query_results": {},
            "query_summary": "",
            "nextflow_config": {},
            "config_reasoning": "",
            "user_decision": "",
            
            # Execute节点输出
            "nextflow_command": "",      # 构建的nextflow命令
            "execution_status": "",      # 执行状态: building/running/completed/failed
            "execution_output": "",      # 执行输出日志
            "execution_result": {},      # 执行结果摘要
        })
