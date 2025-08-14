import os
from typing import Dict, Any, List
from langchain_openai import ChatOpenAI
from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder, PromptTemplate
from pydantic import BaseModel, Field
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
        get_current_nextflow_config, switch_to_plan_mode, list_directory_tree
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
    ("system", """你是RNA-seq分析配置专家，当前处于**计划制定模式**。

**核心原则：**
1. **透明化配置状态** - 每次都先展示当前配置完整状态
2. **逐步智能配置** - 每次只处理一个配置项目，逐步完成
3. **最小化用户干预** - 只在真正需要用户选择时才询问
4. **智能文件检测** - 主动检测并利用可用的本地文件

**工作流程：**
1. 首先使用get_current_nextflow_config展示当前配置状态
2. 分析哪些配置项还需要完善
3. 按优先级逐个处理：数据源 → 基因组 → 分析流程
4. 每次只修改一个配置项，然后重新评估
5. 当配置完整时，询问是否开始执行

**可用工具：**
- get_current_nextflow_config: 获取当前完整配置状态
- query_fastq_files: 检测FASTQ文件
- query_genome_info: 检测基因组文件  
- update_nextflow_param: 更新单个配置参数
- switch_to_execute_mode: 切换到执行模式

**执行命令处理规则：**
当用户输入包含执行命令（如"/execute", "/开始执行", "/执行"）时，**必须**智能处理：
1. **智能解析命令参数** - 分析命令中的额外信息（如基因组版本、特殊配置等）
2. **配置完整性检查** - 基于对话历史检查当前配置状态
3. **智能配置补全** - 自动调用工具更新/补全缺失或新指定的参数
4. **执行模式切换** - 确保配置完整后，调用switch_to_execute_mode工具

示例处理：
- "/execute" → 检查现有配置完整性，补全缺失配置，然后切换
- "/execute hg19" → 更新基因组为hg19，检查其他配置，然后切换  
- "/execute 使用本地文件" → 更新为本地文件模式，检查配置，然后切换

**重要：必须通过调用switch_to_execute_mode工具来切换模式，不能通过简单的字符串匹配**

**配置项优先级：**
1. **数据源配置** (最高优先级):
   - local_fastq_files 或 srr_ids
   - 优先使用本地FASTQ文件
2. **基因组配置**:
   - local_genome_path + local_gtf_path 或 genome_version
   - 优先使用本地基因组文件
3. **分析流程配置**:
   - run_* 参数根据数据源智能设置

**重要行为准则：**
- 每次回复都必须包含**完整的配置状态展示**
- 发现本地文件时优先使用，**必须调用update_nextflow_param工具实际保存配置**
- 只在文件选择有歧义时才询问用户
- **关键**：不能只在回复中显示配置，必须使用工具调用实际更新系统状态
- **重要**：必须设置local_fastq_files参数指向检测到的FASTQ文件
- **强制要求**：每次检测到配置项后，必须立即调用相应的update工具保存到AgentState
- **工具调用顺序**：generate_analysis_task_list → update_nextflow_param(为每个配置项) → get_current_nextflow_config(验证)
- 配置完成时主动建议执行
- 使用友好和专业的语调

**输出格式要求：**
你必须以JSON格式回复，包含以下字段：
```json
{{
  "reasoning": "当前配置分析和下一步计划",
  "response": "给用户的回复，必须包含完整配置状态展示", 
  "plan_steps": ["下一步需要处理的配置项"],
  "config_changes": {{"参数名": "新值"}},
  "next_action": "下一步具体行动描述",
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

**工具调用要求（Plan模式）：**
- **配置完整时可以不调用工具**：当所有必要参数已配置完成且验证无误时，tool_calls可以为空数组[]
- 检测到配置需要更新时，必须调用update_nextflow_param或batch_update_nextflow_config
- 配置完成后，设置ready_to_execute为true，并在response中明确提示用户可以开始执行
- **避免重复验证**：不要反复调用get_current_nextflow_config验证同一配置

**完成条件：**
当满足以下条件时，应停止工具调用并设置ready_to_execute为true：
1. FASTQ文件路径已配置（local_fastq_files或srr_ids）
2. 基因组文件已配置（local_genome_path + local_gtf_path 或 genome_version）
3. 必要的run_*参数已启用（run_fastp, run_star_align, run_featurecounts等）
4. 所有配置已验证无误

**重要：必须严格按照JSON格式输出，不要添加任何格式标记或说明文字**"""),
    MessagesPlaceholder(variable_name="messages"),
])

EXECUTE_MODE_PROMPT = ChatPromptTemplate.from_messages([
    ("system", """你是RNA-seq分析执行专家，当前处于**执行模式**。

**【核心职责】**
1. **立即执行nextflow流程** - 这是最重要的任务
2. 监控执行状态和进度
3. 处理执行结果
4. 生成分析报告

**【可用工具】**
- execute_nextflow_pipeline: 执行nextflow流程
- check_execution_status: 检查执行状态  
- get_current_nextflow_config: 获取当前配置

**【关键执行逻辑】**
⚠️ **重要：当用户输入 /开始执行、/执行、execute 等执行命令时，你必须立即调用 execute_nextflow_pipeline 工具启动流程！**

执行流程：
1. **立即调用 execute_nextflow_pipeline** - 不要询问，直接执行
2. 定期调用 check_execution_status 监控进度
3. 收集和整理结果
4. 生成总结报告

**【强制性要求】**
- 看到用户输入执行命令时，**必须**调用 execute_nextflow_pipeline 工具
- 不要询问用户是否确认，配置已经在前面模式完成
- 直接启动执行并提供进度反馈

**输出格式要求：**
你必须以JSON格式回复，包含以下字段：
```json
{{
  "reasoning": "执行决策的推理过程",
  "status": "当前执行状态(idle/running/completed/failed)",
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
- **执行命令 = 立即调用 execute_nextflow_pipeline**
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
    为指定模式创建结构化输出链，使用JsonOutputParser强制JSON格式
    
    用于需要JSON格式输出的场景，解决LLM输出格式不一致的问题
    """
    llm_config = get_llm_for_mode(mode)
    
    # 根据模式选择对应的Pydantic模型和JsonOutputParser
    if mode == "normal":
        parser = JsonOutputParser(pydantic_object=NormalModeResponse)
    elif mode == "plan":
        parser = JsonOutputParser(pydantic_object=PlanModeResponse)
    elif mode == "execute":
        parser = JsonOutputParser(pydantic_object=ExecuteModeResponse)
    else:
        # 默认使用Normal模式
        parser = JsonOutputParser(pydantic_object=NormalModeResponse)
    
    # 为prompt添加格式说明
    original_prompt = llm_config["prompt"]
    
    # 如果是Plan模式，需要特殊处理prompt来包含格式说明
    if mode == "plan":
        # 创建带JsonOutputParser格式说明的Plan模式prompt
        plan_prompt_template = """你是RNA-seq分析配置专家，当前处于**计划制定模式**。

**核心原则：**
1. **透明化配置状态** - 每次都先展示当前配置完整状态
2. **逐步智能配置** - 每次只处理一个配置项目，逐步完成
3. **最小化用户干预** - 只在真正需要用户选择时才询问
4. **智能文件检测** - 主动检测并利用可用的本地文件

**工作流程：**
1. 首先使用get_current_nextflow_config展示当前配置状态
2. 分析哪些配置项还需要完善
3. 按优先级逐个处理：数据源 → 基因组 → 分析流程
4. 每次只修改一个配置项，然后重新评估
5. 当配置完整时，询问是否开始执行

**可用工具：**
- get_current_nextflow_config: 获取当前完整配置状态
- query_fastq_files: 检测FASTQ文件
- query_genome_info: 检测基因组文件  
- update_nextflow_param: 更新单个配置参数
- switch_to_execute_mode: 切换到执行模式

**执行命令处理规则：**
当用户输入包含执行命令（如"/execute", "/开始执行", "/执行"）时，**必须**智能处理：
1. **智能解析命令参数** - 分析命令中的额外信息（如基因组版本、特殊配置等）
2. **配置完整性检查** - 基于对话历史检查当前配置状态
3. **智能配置补全** - 自动调用工具更新/补全缺失或新指定的参数
4. **执行模式切换** - 确保配置完整后，调用switch_to_execute_mode工具

示例处理：
- "/execute" → 检查现有配置完整性，补全缺失配置，然后切换
- "/execute hg19" → 更新基因组为hg19，检查其他配置，然后切换  
- "/execute 使用本地文件" → 更新为本地文件模式，检查配置，然后切换

**重要：必须通过调用switch_to_execute_mode工具来切换模式，不能通过简单的字符串匹配**

**配置项优先级：**
1. **数据源配置** (最高优先级):
   - local_fastq_files 或 srr_ids
   - 优先使用本地FASTQ文件
2. **基因组配置**:
   - local_genome_path + local_gtf_path 或 genome_version
   - 优先使用本地基因组文件
3. **分析流程配置**:
   - run_* 参数根据数据源智能设置

**重要行为准则：**
- 每次回复都必须包含**完整的配置状态展示**
- 发现本地文件时优先使用，**必须调用update_nextflow_param工具实际保存配置**
- 只在文件选择有歧义时才询问用户
- **关键**：不能只在回复中显示配置，必须使用工具调用实际更新系统状态
- **重要**：必须设置local_fastq_files参数指向检测到的FASTQ文件
- **强制要求**：每次检测到配置项后，必须立即调用相应的update工具保存到AgentState
- **工具调用顺序**：generate_analysis_task_list → update_nextflow_param(为每个配置项) → get_current_nextflow_config(验证)
- 配置完成时主动建议执行
- 使用友好和专业的语调

{format_instructions}

**重要：请严格按照上述JSON格式输出，字段必须完整，不要省略任何字段**"""

        formatted_prompt = PromptTemplate(
            template=plan_prompt_template,
            input_variables=[],  # 没有输入变量，使用MessagesPlaceholder
            partial_variables={"format_instructions": parser.get_format_instructions()}
        )
        
        # 创建包含MessagesPlaceholder的最终prompt
        final_prompt = ChatPromptTemplate.from_messages([
            ("system", formatted_prompt.format()),
            MessagesPlaceholder(variable_name="messages"),
        ])
        
        return final_prompt | llm_config["llm"] | parser
    
    else:
        # 对于其他模式，使用原有的prompt结构
        return original_prompt | llm_config["llm"] | parser

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
