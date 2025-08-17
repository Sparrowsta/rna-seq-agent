import os
from typing import Dict, Any, List
from langchain_deepseek import ChatDeepSeek
from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder, PromptTemplate
import json
from pydantic import BaseModel, Field, field_validator
from .prompt import (
    NORMAL_MODE_PROMPT, PLAN_MODE_PROMPT, EXECUTE_MODE_PROMPT,
    get_prompt_template, get_structured_plan_prompt, PLAN_MODE_STRUCTURED_PROMPT_CONTENT
)
from .tools import (
    query_fastq_files, query_genome_info, add_new_genome,
    update_nextflow_param, batch_update_nextflow_config,
    switch_to_plan_mode, switch_to_execute_mode,
    execute_nextflow_pipeline, check_execution_status, get_current_nextflow_config,
    list_directory_tree, generate_analysis_task_list
)

# ============================================================================
# 工具配置 - 遵循单一职责原则
# ============================================================================

# 所有可用工具
ALL_TOOLS = [
    query_fastq_files, query_genome_info, add_new_genome,
    update_nextflow_param, batch_update_nextflow_config,
    switch_to_plan_mode, switch_to_execute_mode,
    execute_nextflow_pipeline, check_execution_status, get_current_nextflow_config,
    list_directory_tree, generate_analysis_task_list
]

# 模式特定工具映射 - 遵循接口隔离原则
MODE_TOOLS = {
    "normal": [
        query_fastq_files, query_genome_info, add_new_genome,
        get_current_nextflow_config, list_directory_tree
    ],
    "plan": [
        query_fastq_files, query_genome_info, update_nextflow_param,
        batch_update_nextflow_config, get_current_nextflow_config,
        switch_to_execute_mode, generate_analysis_task_list
    ],
    "execute": [
        execute_nextflow_pipeline, check_execution_status,
        get_current_nextflow_config
    ]
}

# ============================================================================
# JSON输出模型定义 - 遵循类型安全原则
# ============================================================================

class ToolCall(BaseModel):
    """工具调用模型"""
    tool_name: str = Field(description="工具名称")
    parameters: Dict[str, Any] = Field(description="工具参数", default_factory=dict)
    reason: str = Field(description="调用原因", default="")

class NormalModeResponse(BaseModel):
    """Normal模式的结构化响应 - 移除tool_calls字段，因为.with_structured_output()内部处理工具调用"""
    reasoning: str = Field(description="分析和推理过程")
    response: str = Field(description="给用户的回复") 
    suggested_actions: List[str] = Field(description="建议的后续操作", default_factory=list)
    need_more_info: bool = Field(description="是否需要更多信息", default=False)

class PlanModeResponse(BaseModel):
    """Plan模式的结构化响应 - 优化为Gemini兼容格式"""
    reasoning: str = Field(
        description="计划制定的推理过程",
        default="正在分析配置需求..."
    )
    plan_steps: List[str] = Field(
        description="分析计划步骤列表", 
        default_factory=list
    )
    config_changes: Dict[str, Any] = Field(
        description="需要修改的nextflow配置参数。必须是有效的JSON对象格式，例如：{{'genome_version': 'hg38', 'run_fastp': true}}。如果没有配置更改，请设置为空对象{{}}。绝对不能使用字符串值如'无'或'已更新'",
        default_factory=dict
    )
    next_action: str = Field(
        description="下一步具体行动建议",
        default="继续配置分析参数"
    )
    ready_to_execute: bool = Field(
        description="是否已准备好执行分析流程", 
        default=False
    )
    tool_calls: List[ToolCall] = Field(
        description="需要调用的工具列表", 
        default_factory=list
    )

class ExecuteModeResponse(BaseModel):
    """Execute模式的结构化响应"""
    reasoning: str = Field(description="执行决策的推理过程")
    status: str = Field(description="当前执行状态")
    progress: str = Field(description="执行进度描述")
    results: Dict[str, Any] = Field(description="执行结果", default_factory=dict)
    next_step: str = Field(description="下一步操作")
    tool_calls: List[ToolCall] = Field(description="需要调用的工具列表", default_factory=list)

# ============================================================================
# LLM配置 - 遵循配置分离原则
# ============================================================================

def create_llm():
    """
    创建DeepSeek LLM实例，使用专门的ChatDeepSeek类
    
    完全替换Gemini，仅使用DeepSeek提供服务
    """
    api_key = os.environ.get("DEEPSEEK_API_KEY")
    
    if not api_key:
        raise ValueError("未找到DEEPSEEK_API_KEY环境变量，请在.env文件中配置")
    
    llm = ChatDeepSeek(
        model="deepseek-chat",
        api_key=api_key,
        temperature=0.1,
    )
    
    return llm

# 基础LLM实例
llm = create_llm()

# ============================================================================
# 模式特定的提示词模板 - 已移至 prompt.py 模块
# ============================================================================

# 注意：所有 prompt 模板已移至 agent/prompt.py 文件
# 这里通过导入方式使用，保持向后兼容

# ============================================================================
# 模式特定的LLM实例 - 遵循策略模式
# ============================================================================

class ModeSpecificLLM:
    """
    模式特定的LLM管理器
    
    遵循策略模式：根据不同模式使用不同的配置
    """
    
    def __init__(self):
        self.base_llm = create_llm()
        self._llm_instances = {}
    
    def get_llm_for_mode(self, mode: str):
        """
        获取特定模式的LLM实例
        
        应用单例模式：每个模式只创建一个LLM实例
        """
        if mode not in self._llm_instances:
            tools = MODE_TOOLS.get(mode, [])
            
            # 使用prompt.py中的get_prompt_template函数统一获取prompt
            prompt_template = get_prompt_template(mode)
            
            self._llm_instances[mode] = {
                "llm": self.base_llm,
                "llm_with_tools": self.base_llm.bind_tools(tools),
                "prompt": prompt_template
            }
        
        return self._llm_instances[mode]

# 全局模式管理器实例
mode_llm_manager = ModeSpecificLLM()

# ============================================================================
# 便捷函数 - 遵循DRY原则
# ============================================================================

def get_llm_for_mode(mode: str):
    """获取指定模式的LLM配置"""
    return mode_llm_manager.get_llm_for_mode(mode)

def create_chain_for_mode(mode: str):
    """
    为指定模式创建处理链
    
    应用链式调用模式：prompt + llm
    """
    llm_config = get_llm_for_mode(mode)
    return llm_config["prompt"] | llm_config["llm_with_tools"]

def create_tool_calling_chain(mode: str):
    """
    创建纯工具调用链 - 双LLM架构第一阶段
    
    专门负责工具调用，不处理结构化输出
    """
    llm_config = get_llm_for_mode(mode)
    return llm_config["prompt"] | llm_config["llm_with_tools"]

def create_structured_response_chain(mode: str):
    """
    创建纯结构化输出链 - 双LLM架构第二阶段
    
    专门负责将工具调用结果格式化为结构化输出
    """
    from langchain_core.prompts import ChatPromptTemplate
    
    # 根据模式选择对应的Pydantic模型
    if mode == "normal":
        pydantic_model = NormalModeResponse
    elif mode == "plan":
        pydantic_model = PlanModeResponse
    elif mode == "execute":
        pydantic_model = ExecuteModeResponse
    else:
        pydantic_model = NormalModeResponse
    
    # 创建专门的结构化输出提示词
    structured_prompt = ChatPromptTemplate.from_messages([
        ("system", f"""你是RNA-seq分析助手。根据用户的原始请求和工具执行结果，生成结构化的响应。

你需要分析工具调用的结果，并将其格式化为用户友好的响应格式。

重要：你的响应必须是有效的JSON格式，包含以下字段：
- reasoning: 分析和推理过程
- response: 给用户的主要回复
- suggested_actions: 建议的后续操作列表
- need_more_info: 是否需要更多信息

确保响应内容准确、专业且有帮助。"""),
        ("user", "原始用户请求：{original_input}\n\n工具执行结果：{tool_results}\n\n请生成结构化响应：")
    ])
    
    # 获取基础LLM（不绑定工具）
    llm_config = get_llm_for_mode(mode)
    base_llm = llm_config["llm"]
    
    # 创建结构化输出LLM
    structured_llm = base_llm.with_structured_output(pydantic_model, method="json_mode")
    
    return structured_prompt | structured_llm

def create_dual_llm_chain_for_mode(mode: str):
    """
    创建双LLM链 - 组合工具调用和结构化输出
    
    第一阶段：工具调用链处理工具执行
    第二阶段：结构化输出链格式化响应
    """
    tool_chain = create_tool_calling_chain(mode)
    response_chain = create_structured_response_chain(mode)
    
    return {
        "tool_chain": tool_chain,
        "response_chain": response_chain
    }

def create_structured_chain_for_mode(mode: str):
    """
    向后兼容函数 - 返回双LLM链配置
    
    保持现有接口，内部使用双LLM架构
    """
    return create_dual_llm_chain_for_mode(mode)

# ============================================================================
# 向后兼容 - 保持现有接口
# ============================================================================

# 为了保持向后兼容，保留原有的变量名
tools = ALL_TOOLS
llm_with_tools = llm.bind_tools(tools)

# 默认prompt（normal模式）
prompt = get_prompt_template("normal")

# ============================================================================
# 配置验证 - 确保系统健壮性
# ============================================================================

def validate_environment():
    """
    验证环境配置
    
    应用KISS原则：简单的环境检查
    """
    required_env_vars = ["OPENAI_API_KEY"]
    missing_vars = []
    
    for var in required_env_vars:
        if not os.environ.get(var):
            missing_vars.append(var)
    
    if missing_vars:
        raise ValueError(f"缺少必需的环境变量: {', '.join(missing_vars)}")
    
    return True

# 在模块加载时验证环境
try:
    validate_environment()
except ValueError as e:
    print(f"警告：{e}")
