import json
from typing import Dict, Any
from ..state import AgentState, PrepareResponse
from ..core import get_shared_llm

def create_prepare_agent():
    """创建Prepare节点的智能配置Agent"""
    llm = get_shared_llm()
    structured_llm = llm.with_structured_output(PrepareResponse, method="json_mode")
    return structured_llm

async def prepare_node(state: AgentState) -> Dict[str, Any]:
    """准备节点 - 优先基于Normal模式传来的用户需求生成配置参数"""
    print(f"⚙️ 开始智能配置分析...")
    
    # 获取所有必要信息
    detection_results = state.query_results or {}
    current_config = state.nextflow_config or {}
    initial_requirements = state.user_requirements or {}
    replan_requirements = state.replan_requirements or {}
    
    print(f"📝 初始配置需求: {initial_requirements}")
    if replan_requirements:
        print(f"🔄 重新规划需求: {replan_requirements}")
    print(f"📊 当前配置状态: {current_config}")
    
    if not detection_results:
        print("⚠️ 未检测到任何数据，无法生成配置")
        return {
            "nextflow_config": current_config,
            "config_reasoning": "未获取到检测数据，保持现有配置",
            "response": "⚠️ 缺少检测数据，无法进行智能配置分析",
            "status": "error"
        }
    
    # 使用LLM综合分析用户需求和检测数据
    prepare_agent = create_prepare_agent()
    
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
3. **填充缺失配置** - 对用户未指定的字段使用系统推荐值

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
注意：不是字典格式，是数组格式！

**必需配置字段：**
- genome_version, species: 基因组相关（优先使用用户指定值）
- qc_tool, align_tool, quant_tool: 工具链（小写，优先使用用户指定值）
- local_fastq_files: 文件路径列表
- paired_end: 是否包含双端数据
- sample_groups: 详细样本配对信息
- run_build_star_index: 索引构建控制
- run_download_genome: 基因组下载控制

**决策说明要求：**
在config_reasoning中以文本格式详细说明：
1. 用户需求如何被直接应用 (初始需求: [initial_requirements], 重新规划需求: [replan_requirements])  
2. 系统检测结果在哪些字段被使用
3. 每个关键配置的最终决策理由

**返回JSON格式字段：**
- nextflow_config: 完整的Nextflow配置参数字典
- config_reasoning: 配置决策理由的详细文本说明（字符串格式，不是嵌套字典）

**config_reasoning格式示例：**
"基于用户需求分析：无特殊要求，采用系统推荐配置。基因组选择：hg19_test因文件完整性最佳。工具选择：fastp+STAR+featureCounts基于可用性检测。FASTQ配对：检测到3个双端样本，生成数组格式sample_groups。索引策略：现有STAR索引完整，无需重建。"

=== 📊 系统检测数据 ===
{json.dumps(detection_results, indent=2, ensure_ascii=False)}

=== ⚙️ 当前配置状态 ===
{json.dumps(current_config, indent=2, ensure_ascii=False)}

**重要提醒：**
1. 重新规划需求优先级高于初始需求，如有冲突以重新规划需求为准
2. 基于检测数据进行FASTQ文件智能配对
3. 对用户未指定的字段使用系统检测推荐值
4. 在reasoning中详细说明用户需求的应用情况"""
        
        # LLM直接输出PrepareResponse格式
        analysis_result = await prepare_agent.ainvoke([{"role": "user", "content": unified_prompt}])
        
        # 检查LLM响应
        if not analysis_result:
            raise Exception("LLM返回空响应")
        
        # 提取结果
        config_params = analysis_result.nextflow_config or {}
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
            "config_reasoning": reasoning,
            "response": f"智能配置分析完成{user_satisfaction_note}\n\n💡 {reasoning}\n\n🔧 生成了 {len(config_params)} 个配置参数",
            "status": "confirm"
        }
        
    except Exception as e:
        print(f"❌ LLM分析失败: {str(e)}")
        return {
            "nextflow_config": current_config,
            "config_reasoning": f"LLM分析失败: {str(e)}",
            "response": f"❌ 配置生成失败: {str(e)}",
            "status": "error"
        }