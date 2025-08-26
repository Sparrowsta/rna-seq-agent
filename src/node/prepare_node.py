import json
from typing import Dict, Any
from ..state import AgentState, PrepareResponse
from ..core import get_shared_llm
from langgraph.prebuilt import create_react_agent

def create_prepare_agent():
    """创建Prepare节点的智能配置Agent"""
    llm = get_shared_llm()
    
    # 使用create_react_agent但不提供tools，纯推理模式
    agent = create_react_agent(
        model=llm,
        tools=[],  # 空工具列表，纯推理
        prompt="你是RNA-seq分析配置专家。请基于用户需求和系统检测数据，生成最优化的Nextflow配置参数。",
        response_format=PrepareResponse
    )
    return agent

async def prepare_node(state: AgentState) -> Dict[str, Any]:
    """准备节点 - 优先基于Normal模式传来的用户需求生成配置参数"""
    print(f"⚙️ 开始智能配置分析...")
    
    # 获取所有必要信息
    detection_results = state.query_results or {}
    current_config = state.nextflow_config or {}
    initial_requirements = state.user_requirements or {}
    replan_requirements = state.replan_requirements or {}
    
    if not detection_results:
        return {
            "nextflow_config": current_config,
            "resource_config": {},
            "config_reasoning": "未获取到检测数据，保持现有配置",
            "response": "⚠️ 缺少检测数据，无法进行智能配置分析",
            "status": "error"
        }
    
    # 使用LLM综合分析用户需求和检测数据
    
    try:
        print("🧠 LLM正在基于用户需求分析检测数据并生成配置...")
        
        # 构建统一的完整prompt
        requirements_section = ""
        if initial_requirements:
            requirements_section += f"\n**初始配置需求: {initial_requirements}**"
        if replan_requirements:
            requirements_section += f"\n**重新规划需求: {replan_requirements}** (优先级更高)"
        
        unified_prompt = f"""你是RNA-seq分析配置专家。请基于用户需求和检测数据生成最优配置。

请以JSON格式返回分析结果。

{requirements_section}

**需求处理策略：**
- 如果同时存在初始需求和重新规划需求，优先满足重新规划需求
- 重新规划需求可以覆盖初始需求中的任何配置项
- 综合考虑所有需求，确保生成的配置满足用户的最终意图

**配置决策优先级：**
1. **重新规划需求绝对优先** - 如存在重新规划需求，优先采用
2. **初始需求作为基础** - 初始需求作为基础配置参考
3. **技术可行性适配** - 确保配置在技术上可行
4. **系统智能推荐** - 仅在用户未指定的配置项使用检测推荐值

**核心任务：**
1. **应用用户配置** - 优先级：重新规划需求 > 初始需求 > 系统推荐
2. **FASTQ文件配对分析** - 基于fastq_analysis进行智能文件配对
3. **资源智能分配** - 基于系统检测和样本规模进行CPU/内存优化
**基因组配置 - 对用户想要使用的基因组进行配置，没有则按照系统推荐。根据基因组是否存在，基因组索引是否构建来调整对应的配置字段：

**关键配置决策逻辑：**
- **run_download_genome**: 
  - 如果基因组文件(FASTA+GTF)都已存在 → 设为 false
  - 如果基因组文件缺失或不完整 → 设为 true
- **run_build_star_index**: 
  - 如果STAR索引目录已存在且包含完整索引文件 → 设为 false  
  - 如果STAR索引不存在或不完整 → 设为 true

**基因组状态检查重点：**
从系统检测数据的genome_analysis中查看：
- 每个基因组的fasta_file.exists和gtf_file.exists状态
- 每个基因组的star_index.exists状态和file_count
- 优先选择文件完整且已有索引的基因组减少处理时间

**示例决策：**
- 如果hg19的FASTA/GTF存在且STAR索引已构建(file_count > 10) → run_download_genome: false, run_build_star_index: false
- 如果用户指定基因组但文件不存在 → run_download_genome: true, run_build_star_index: true

**FASTQ配对分析：**
从fastq_analysis.file_paths分析文件名模式并使用完整路径：
- 双端：sample_1.fastq.gz + sample_2.fastq.gz  
- 双端：sample_R1.fastq.gz + sample_R2.fastq.gz
- 单端：sample.fastq.gz
- **重要：必须使用file_paths中的完整路径，如"fastq/SRR17469061_1.fastq.gz"**

**sample_groups格式要求（重要）：**
必须生成数组格式，每个元素包含sample_id、read1、read2字段，使用完整文件路径：
[
  {{"sample_id": "SRR17469061", "read1": "fastq/SRR17469061_1.fastq.gz", "read2": "fastq/SRR17469061_2.fastq.gz"}},
  {{"sample_id": "SRR17469059", "read1": "fastq/SRR17469059_1.fastq.gz", "read2": "fastq/SRR17469059_2.fastq.gz"}}
]
注意：不是字典格式，是数组格式

**必需配置字段：**
- genome_version, species: 基因组相关（优先使用用户指定值）
- qc_tool, align_tool, quant_tool: 工具链（小写，优先使用用户指定值）
- paired_end: 是否包含双端数据
- sample_groups: 详细样本配对信息
- run_build_star_index: 索引构建控制
- run_download_genome: 基因组下载控制

**资源配置决策（核心新增功能）：**
基于system_resources检测数据进行智能资源分配：

**资源分配策略：**
1. **系统资源评估** - 从system_resources获取：
   - total_cpus: 总CPU核心数
   - total_memory_gb: 总内存（GB）
   - disk_space_gb: 磁盘空间（GB）

2. **样本规模分析** - 从fastq_analysis获取：
   - total_files_found: FASTQ文件总数
   - file_size_summary: 文件大小信息
   - 推测数据处理复杂度

3. **智能资源分配规则：**
   - **CPU分配原则**: 不超过total_cpus的80%，预留20%给系统
   - **内存分配原则**: 基于文件大小和进程类型智能调整
   - **进程优先级**: prepare_star_index > run_alignment > run_quality_control > run_quantification
   - **工具基本要求**: 如果使用“star”工具，则内存分配至少要32GB
4. **具体分配策略：**
   ```
   prepare_star_index: max(总CPU*0.6, 4核), 至少32GB内存(STAR工具要求)
   run_alignment: max(总CPU*0.5, 4核), 至少32GB内存(STAR工具要求)
   run_quality_control: max(总CPU*0.4, 2核), max(总内存*0.2, 8GB)
   run_quantification: max(总CPU*0.3, 2核), max(总内存*0.15, 6GB)
   download进程: 2核, 4GB (IO密集型，CPU需求低)
   ```

5. **STAR内存要求特别处理：**
   - prepare_star_index和run_alignment进程强制最少32GB内存
   - 如果系统总内存<32GB，仍设置为32GB但会警告用户内存不足风险
   - 内存不足时建议：增加虚拟内存、减少其他程序、使用更小的基因组

6. **文件大小适配：**
   - 大文件(>2GB): 内存需求 × 1.5倍
   - 小文件(<500MB): 内存需求 × 0.8倍
   - 样本数量>5个: CPU需求 × 1.2倍（并行处理优势）

**resource_config输出格式：**
必须生成字典格式，包含每个进程的资源配置：
```json
{{
  "prepare_star_index": {{"cpus": 8, "memory": "32 GB", "reasoning": "STAR索引构建强制32GB，系统内存充足"}},
  "run_alignment": {{"cpus": 6, "memory": "32 GB", "reasoning": "STAR比对强制32GB，系统内存充足"}}, 
  "run_quality_control": {{"cpus": 4, "memory": "12 GB", "reasoning": "质控处理中等资源需求"}},
  "run_quantification": {{"cpus": 3, "memory": "8 GB", "reasoning": "定量分析轻量级处理"}},
  "download_genome_fasta": {{"cpus": 2, "memory": "4 GB", "reasoning": "IO密集型下载任务"}},
  "download_genome_gtf": {{"cpus": 2, "memory": "4 GB", "reasoning": "IO密集型下载任务"}}
}}
```

**决策说明要求：**
在config_reasoning中以文本格式详细说明：
1. 用户需求如何被直接应用 (初始需求: [initial_requirements], 重新规划需求: [replan_requirements])  
2. 基因组索引决策的详细分析 - **必须明确说明 run_download_genome 和 run_build_star_index 的设置理由**
3. **资源分配决策过程** - 系统资源检测结果、样本规模评估、资源分配策略应用
4. 系统检测结果在哪些字段被使用
5. 每个关键配置的最终决策理由

**返回JSON格式字段：**
- nextflow_config: 完整的Nextflow配置参数字典
- resource_config: 各进程的CPU和内存资源配置字典
- config_reasoning: 配置决策理由的详细文本说明（字符串格式，不是嵌套字典）

**config_reasoning格式示例：**
"基于用户需求分析：用户指定hg19基因组进行分析。基因组配置检查：hg19 FASTA/GTF文件已存在，STAR索引已构建(15个文件)，因此设置 run_download_genome: false, run_build_star_index: false。工具选择：fastp+STAR+featureCounts基于用户要求和系统可用性。FASTQ配对：检测到3个双端样本，生成数组格式sample_groups。最终配置优化了处理时间，跳过不必要的下载和索引构建。"

=== 📊 系统检测数据 ===
{json.dumps(detection_results, indent=2, ensure_ascii=False)}

=== ⚙️ 当前配置状态 ===
{json.dumps(current_config, indent=2, ensure_ascii=False)}

**重要提醒：**
1. 重新规划需求优先级高于初始需求，如有冲突以重新规划需求为准
2. 基于检测数据进行FASTQ文件智能配对
3. 对用户未指定的字段使用系统检测推荐值
4. 在reasoning中详细说明用户需求的应用情况"""
        
        # 使用create_react_agent调用方式
        agent_executor = create_prepare_agent()
        messages_input = {"messages": [{"role": "user", "content": unified_prompt}]}
        
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")
        
        if structured_response:
            analysis_result = structured_response
        else:
            raise Exception("Agent未返回预期的结构化响应")
        
        # 检查LLM响应
        if not analysis_result:
            raise Exception("LLM返回空响应")
        
        # 提取结果
        config_params = analysis_result.nextflow_config or {}
        resource_params = analysis_result.resource_config or {}
        reasoning = analysis_result.config_reasoning or "基于用户需求和检测数据的智能分析"
        
        print(f"✅ 配置生成完成，严格遵循用户需求")
        
        # 合并配置参数（新配置优先）
        final_config = current_config.copy()
        final_config.update(config_params)
        
        # 构建需求满足情况说明
        user_satisfaction_note = ""
        if initial_requirements or replan_requirements:
            satisfaction_parts = []
            if initial_requirements:
                satisfaction_parts.append(f"初始需求: {initial_requirements}")
            if replan_requirements:
                satisfaction_parts.append(f"重新规划需求: {replan_requirements} (已优先应用)")
            user_satisfaction_note = f"\n\n🎯 **用户需求满足情况：**\n" + "\n".join(satisfaction_parts)
        
        return {
            "nextflow_config": final_config,
            "resource_config": resource_params,
            "config_reasoning": reasoning,
            "response": f"智能配置分析完成{user_satisfaction_note}\n\n💡 {reasoning}\n\n🔧 生成了 {len(config_params)} 个配置参数\n🖥️ 生成了 {len(resource_params)} 个进程的资源配置",
            "status": "confirm"
        }
        
    except Exception as e:
        print(f"❌ LLM分析失败: {str(e)}")
        return {
            "nextflow_config": current_config,
            "resource_config": {},
            "config_reasoning": f"LLM分析失败: {str(e)}",
            "response": f"❌ 配置生成失败: {str(e)}",
            "status": "error"
        }