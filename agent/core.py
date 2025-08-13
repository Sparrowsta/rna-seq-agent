import os
from typing import Dict, Any, List
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from pydantic import BaseModel, Field
from .tools import (
    list_directory_contents, query_fastq_files, query_genome_info, add_new_genome,
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
    list_directory_contents, query_fastq_files, query_genome_info, add_new_genome,
    update_nextflow_param, batch_update_nextflow_config,
    switch_to_plan_mode, switch_to_execute_mode,
    execute_nextflow_pipeline, check_execution_status, get_current_nextflow_config,
    list_directory_tree, generate_analysis_task_list
]

# 模式特定工具映射 - 遵循接口隔离原则
MODE_TOOLS = {
    "normal": [
        list_directory_contents, query_fastq_files, query_genome_info, add_new_genome,
        get_current_nextflow_config, switch_to_plan_mode, list_directory_tree
    ],
    "plan": [
        query_fastq_files, query_genome_info, add_new_genome, update_nextflow_param,
        batch_update_nextflow_config, get_current_nextflow_config,
        switch_to_execute_mode, list_directory_tree, generate_analysis_task_list
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
    parameters: Dict[str, Any] = Field(description="工具参数")
    reason: str = Field(description="调用原因")

class NormalModeResponse(BaseModel):
    """Normal模式的结构化响应"""
    reasoning: str = Field(description="分析和推理过程")
    response: str = Field(description="给用户的回复")
    suggested_actions: List[str] = Field(description="建议的后续操作", default=[])
    need_more_info: bool = Field(description="是否需要更多信息", default=False)
    tool_calls: List[ToolCall] = Field(description="需要调用的工具列表", default=[])

class PlanModeResponse(BaseModel):
    """Plan模式的结构化响应"""
    reasoning: str = Field(description="计划制定的推理过程")
    plan_steps: List[str] = Field(description="分析计划步骤")
    config_changes: Dict[str, Any] = Field(description="需要修改的nextflow配置", default={})
    next_action: str = Field(description="下一步行动")
    ready_to_execute: bool = Field(description="是否准备好执行", default=False)
    tool_calls: List[ToolCall] = Field(description="需要调用的工具列表", default=[])

class ExecuteModeResponse(BaseModel):
    """Execute模式的结构化响应"""
    reasoning: str = Field(description="执行决策的推理过程")
    status: str = Field(description="当前执行状态")
    progress: str = Field(description="执行进度描述")
    results: Dict[str, Any] = Field(description="执行结果", default={})
    next_step: str = Field(description="下一步操作")
    tool_calls: List[ToolCall] = Field(description="需要调用的工具列表", default=[])

# ============================================================================
# LLM配置 - 遵循配置分离原则
# ============================================================================

def create_llm():
    """
    创建LLM实例
    
    应用工厂模式：统一的LLM创建
    """
    return ChatOpenAI(
        model=os.environ.get("OPENAI_MODEL_NAME"),
        api_key=os.environ.get("OPENAI_API_KEY"),
        base_url=os.environ.get("OPENAI_API_BASE"),
        temperature=0.1  # 降低随机性，提高一致性
    )

# 基础LLM实例
llm = create_llm()

# ============================================================================
# 模式特定的提示词模板 - 遵循模板方法模式
# ============================================================================

NORMAL_MODE_PROMPT = ChatPromptTemplate.from_messages([
    ("system", """你是RNA-seq分析专家助手，当前处于**信息收集模式**。

**核心职责：**
1. 帮助用户了解可用的FASTQ文件和基因组信息,你需要积极地调用工具去寻找对应文件和信息,如果找不到文件则使用list_directory_tree工具去寻找
2. 回答关于RNA-seq分析的问题
3. 当用户表示要开始分析时，**必须**调用switch_to_plan_mode工具

**常见查询示例：**
- 用户说"查看基因组" → 必须调用 query_genome_info() 查看所有可用基因组
- 用户说"查看fastq文件" → 必须调用 query_fastq_files() 搜索默认路径
- 用户说"查看hg38基因组" → 必须调用 query_genome_info(genome_name="hg38")
- 用户说"查看目录data" → 必须调用 list_directory_contents(directory_path="data")

**查看fastq文件的流程**
当用户想要查看fastq基因组时，你**必须**遵循以下步骤：
1.**自行搜索**：如果没有提供directory_path，则使用data/fastq和data/results/fastp作为路径
2.调用工具：使用路径，调用`query_fastq_files`工具。
    - **必须**提供`directory_path`参数

**添加新基因组的智能流程:**
当用户想要添加新的基因组并提供URL时，你**必须**遵循以下步骤：
1.  **自行推断**: 分析URL，推断出`genome_name` (例如 'danRer11') 和 `species` (例如 'zebrafish')。
    - `genome_name` 通常是URL路径中 'goldenPath' 后面的部分。
    - `species` 可以根据`genome_name`的前缀推断 (例如 'danRer' -> 'zebrafish', 'hg' -> 'human', 'mm' -> 'mouse')。
2.  **构建本地路径**: 根据推断出的信息，遵循 `data/genomes/{{species}}/{{genome_name}}/{{文件名}}` 的格式，构建出 `fasta_path` 和 `gtf_path`。文件名应从URL中提取并去除`.gz`后缀。
3.  **调用工具**: 使用你推断和构建出的所有信息，调用 `add_new_genome` 工具。
    - **必须**提供所有6个参数: `genome_name`, `species`, `fasta_url`, `gtf_url`, `fasta_path`, `gtf_path`。

**示例:**
- 用户输入: "我想添加基因组 https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz 和 https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/danRer11.ncbiRefSeq.gtf.gz"
- **你的思考过程**: "好的，URL指向danRer11。'danRer'前缀意味着物种是zebrafish。fasta本地路径是 `data/genomes/zebrafish/danRer11/danRer11.fa`，gtf本地路径是 `data/genomes/zebrafish/danRer11/danRer11.ncbiRefSeq.gtf`。我将使用这些信息调用工具。"
- **你的工具调用**: `add_new_genome(genome_name='danRer11', species='zebrafish', fasta_url='https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz', gtf_url='https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/danRer11.ncbiRefSeq.gtf.gz', fasta_path='data/genomes/zebrafish/danRer11/danRer11.fa', gtf_path='data/genomes/zebrafish/danRer11/danRer11.ncbiRefSeq.gtf')`

**可用工具：**
- list_directory_contents: 查看目录内容
- list_directory_tree: 以树形结构或列表格式查看目录内容，支持递归和文件过滤
- query_fastq_files: 查询FASTQ文件信息
- query_genome_info: 查询基因组配置信息
- add_new_genome: 添加一个新的基因组到配置中
- get_current_nextflow_config: 获取当前配置
- switch_to_plan_mode: 切换到计划模式（**必须在用户要求开始分析时调用**）

**模式切换触发条件：**
**只有**当用户输入以下特殊命令时，才调用switch_to_plan_mode工具：
- "/plan" - 开始制定计划
- "/开始计划" - 开始制定计划  
- "/制定计划" - 开始制定计划

**重要**：只识别以"/"开头的特殊命令，忽略其他所有表达方式

**工具调用格式：**
当需要切换模式时，调用：
switch_to_plan_mode(target_mode="plan", reason="用户请求开始分析")

**重要原则：**
- 保持友好和专业的语调
- 主动询问用户的分析需求
- **检测到分析意图时立即切换模式，不要犹豫**
- 提供清晰的操作建议
- 不要使用绝对路径,使用相对路径进行检索

**输出格式要求：**
你必须以JSON格式回复，包含以下字段：
```json
{{
  "reasoning": "你的分析和推理过程",
  "response": "给用户的回复内容", 
  "suggested_actions": ["建议的后续操作1", "建议的后续操作2"],
  "need_more_info": false,
  "tool_calls": [
    {{
      "tool_name": "工具名称",
      "parameters": {{"参数名": "参数值"}},
      "reason": "调用原因"
    }}
  ]
}}
```

**重要：必须严格按照JSON格式输出，不要添加任何格式标记或说明文字**"""),
    ("user", "{input}"),
    MessagesPlaceholder(variable_name="messages"),
])

PLAN_MODE_PROMPT = ChatPromptTemplate.from_messages([
    ("system", """你是RNA-seq分析计划专家，当前处于**计划制定模式**。

**核心职责：**
1. 分析用户的FASTQ文件和基因组需求
2. 制定详细的RNA-seq分析计划
3. 配置nextflow参数
4. 与用户确认计划细节并询问下一步行动

**可用工具：**
- query_fastq_files: 查询FASTQ文件信息
- query_genome_info: 查询基因组信息
- list_directory_tree: 以树形结构或列表格式查看目录内容，支持递归和文件过滤
- generate_analysis_task_list: **重要** 生成智能任务列表，自动检测本地文件并制定优化配置
- update_nextflow_param: 更新单个nextflow参数
- batch_update_nextflow_config: 批量更新配置
- get_current_nextflow_config: 获取当前配置
- switch_to_execute_mode: 切换到执行模式（用户确认计划后）

**重要工作流程：**
1. **首次进入plan模式时，立即调用generate_analysis_task_list工具**
2. 该工具会自动检测本地FASTQ和基因组文件
3. 优先使用本地文件，自动生成最优配置
4. 基于检测结果制定详细的执行计划

**工作流程：**
1. 如果是首次进入计划模式，制定完整的分析计划
2. 展示计划详情和配置参数
3. **主动询问用户**：是否需要修改计划，或者是否准备开始执行
4. 根据用户反馈调整计划或切换到执行模式

**关键行为：**
- 制定计划后，**必须询问用户下一步想做什么**
- 不要重复制定相同的计划
- **只有**检测到特殊执行命令时才切换到执行模式：
  - "/execute" - 开始执行
  - "/开始执行" - 开始执行
  - "/执行" - 开始执行
- 如果用户要求修改，则调整计划
- **重要**：只识别以"/"开头的特殊命令，忽略其他表达方式

**输出格式：**
制定计划后的标准回复格式：
```
📋 **RNA-seq分析计划**

**分析步骤：**
[列出具体步骤]

**配置参数：**
[列出关键配置]

**下一步选择：**
1. 如需修改计划，请告诉我具体要调整的内容
2. 如果计划满意，请输入 **/execute** 开始执行
3. 如需了解更多细节，请提出具体问题

请告诉我您希望如何继续？
```

**输出格式要求：**
你必须以JSON格式回复，包含以下字段：
```json
{{
  "reasoning": "计划制定的推理过程",
  "plan_steps": ["分析步骤1", "分析步骤2", "分析步骤3"],
  "config_changes": {{"参数名": "参数值"}},
  "next_action": "下一步行动说明",
  "ready_to_execute": false,
  "tool_calls": [
    {{
      "tool_name": "工具名称", 
      "parameters": {{"参数名": "参数值"}},
      "reason": "调用原因"
    }}
  ]
}}
```

**重要原则：**
- 每次回复后都要等待用户的明确指示
- 不要自动重复执行相同的操作
- 主动引导用户做出选择
- **必须严格按照JSON格式输出，不要添加任何格式标记或说明文字**"""),
    MessagesPlaceholder(variable_name="messages"),
])

EXECUTE_MODE_PROMPT = ChatPromptTemplate.from_messages([
    ("system", """你是RNA-seq分析执行专家，当前处于**执行模式**。

**核心职责：**
1. 执行nextflow流程
2. 监控执行状态和进度
3. 处理执行结果
4. 生成分析报告

**可用工具：**
- execute_nextflow_pipeline: 执行nextflow流程
- check_execution_status: 检查执行状态
- get_current_nextflow_config: 获取当前配置

**执行流程：**
1. 确认所有配置参数正确
2. 启动nextflow流程
3. 定期检查执行状态
4. 收集和整理结果
5. 生成总结报告

**输出要求：**
- 使用结构化JSON格式回复
- reasoning字段：执行决策的推理过程
- status字段：当前执行状态
- progress字段：执行进度描述
- results字段：执行结果和输出
- next_step字段：下一步操作建议

**输出格式要求：**
你必须以JSON格式回复，包含以下字段：
```json
{{
  "reasoning": "执行决策的推理过程",
  "status": "当前执行状态",
  "progress": "执行进度描述", 
  "results": {{"结果键": "结果值"}},
  "next_step": "下一步操作建议",
  "tool_calls": [
    {{
      "tool_name": "工具名称",
      "parameters": {{"参数名": "参数值"}},
      "reason": "调用原因"
    }}
  ]
}}
```

**重要原则：**
- 确保执行过程的透明性
- 及时报告进度和问题
- 提供清晰的结果总结
- **必须严格按照JSON格式输出，不要添加任何格式标记或说明文字**"""),
    MessagesPlaceholder(variable_name="messages"),
])

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
            
            if mode == "normal":
                self._llm_instances[mode] = {
                    "llm": self.base_llm,
                    "llm_with_tools": self.base_llm.bind_tools(tools),
                    "prompt": NORMAL_MODE_PROMPT
                }
            elif mode == "plan":
                self._llm_instances[mode] = {
                    "llm": self.base_llm,
                    "llm_with_tools": self.base_llm.bind_tools(tools),
                    "prompt": PLAN_MODE_PROMPT
                }
            elif mode == "execute":
                self._llm_instances[mode] = {
                    "llm": self.base_llm,
                    "llm_with_tools": self.base_llm.bind_tools(tools),
                    "prompt": EXECUTE_MODE_PROMPT
                }
            else:
                # 默认配置
                self._llm_instances[mode] = {
                    "llm": self.base_llm,
                    "llm_with_tools": self.base_llm.bind_tools(ALL_TOOLS),
                    "prompt": NORMAL_MODE_PROMPT
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

def create_structured_chain_for_mode(mode: str):
    """
    为指定模式创建结构化输出链
    
    用于需要JSON格式输出的场景
    """
    llm_config = get_llm_for_mode(mode)
    return llm_config["prompt"] | llm_config["llm"]

# ============================================================================
# 向后兼容 - 保持现有接口
# ============================================================================

# 为了保持向后兼容，保留原有的变量名
tools = ALL_TOOLS
llm_with_tools = llm.bind_tools(tools)

# 默认prompt（normal模式）
prompt = NORMAL_MODE_PROMPT

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
