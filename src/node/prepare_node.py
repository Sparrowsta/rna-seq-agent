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
    """准备节点 - 综合用户需求和检测数据生成配置参数"""
    print(f"⚙️ 开始智能配置分析...")
    
    # 获取所有必要信息
    detection_results = state.query_results or {}
    current_config = state.nextflow_config or {}
    user_requirements = state.user_requirements or ""
    
    if user_requirements:
        print(f"📝 用户需求: {user_requirements}")
    
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
        print("🧠 LLM正在分析检测数据并生成配置...")
        
        # 构建综合分析的系统消息
        system_message = """你是RNA-seq分析配置专家。综合用户需求和检测数据生成最优配置。

**核心任务：智能FASTQ文件配对分析**
从fastq_analysis.file_paths中的文件列表，根据文件名模式智能识别：
1. 样本分组：提取样本ID (移除_1/_2/_R1/_R2等后缀)
2. 配对关系：判断单端(single-end)还是双端(paired-end)测序
3. 生成sample_groups结构：为每个样本指定read1/read2文件

常见文件名模式：
- 双端：sample_1.fastq.gz + sample_2.fastq.gz
- 双端：sample_R1.fastq.gz + sample_R2.fastq.gz  
- 单端：sample.fastq.gz (没有配对文件)

配置决策优先级：
1. **用户明确需求优先** - 如用户指定基因组版本，必须按要求设置
2. **技术可行性考虑** - 如果用户需求的资源不存在，说明需要下载
3. **系统推荐默认值** - 在用户没有明确要求时使用检测到的可用资源

重要配置字段：
- genome_version, species: 基因组相关配置
- qc_tool, align_tool, quant_tool: 工具选择（必须使用小写：fastp, star, featurecounts）  
- local_fastq_files: 原始文件路径列表
- paired_end: 整体是否包含双端测序
- sample_groups: 每个样本的详细配对信息
- run_build_star_index：当STAR索引没有建立的时候，要启动本地构建基因组
- run_download_genome：当所需的基因组没有在本地时，需要通过url下载

sample_groups格式示例：
[
  {
    "sample_id": "SRR17469061",
    "read1": "fastq/SRR17469061_1.fastq.gz",
    "read2": "fastq/SRR17469061_2.fastq.gz", 
    "is_paired": true
  },
  {
    "sample_id": "sample_single",
    "read1": "fastq/sample_single.fastq.gz",
    "read2": null,
    "is_paired": false
  }
]

返回JSON格式：
{
  "nextflow_config": {
    "genome_version": "hg38",
    "local_fastq_files": ["file1.fastq.gz", "file2.fastq.gz"],
    "paired_end": true,
    "sample_groups": [样本配对数组]
  },
  "config_reasoning": "详细分析说明"
}"""
        
        # 构建包含所有信息的用户消息
        user_message_parts = [
            "请综合以下信息生成最优配置：",
            "",
            "=== 检测数据 ===",
            json.dumps(detection_results, indent=2, ensure_ascii=False),
            "",
            f"=== 当前配置 ===", 
            json.dumps(current_config, indent=2, ensure_ascii=False)
        ]
        
        if user_requirements:
            user_message_parts.extend([
                "",
                f"=== 用户明确需求 ===",
                f"**{user_requirements}**",
                "",
                "注意：用户需求应优先满足，即使检测数据显示相关资源不存在，也要按用户要求配置，并在reasoning中说明解决方案。"
            ])
        
        user_message_parts.extend([
            "",
            "请分析所有信息，生成既满足用户需求又考虑技术可行性的配置参数。"
        ])
        
        user_message = "\n".join(user_message_parts)
        
        messages = [
            {"role": "system", "content": system_message},
            {"role": "user", "content": user_message}
        ]
        
        # LLM直接输出PrepareResponse格式
        analysis_result = prepare_agent.invoke(messages)
        
        # 检查LLM响应
        if not analysis_result:
            raise Exception("LLM返回空响应")
        
        # 提取结果
        config_params = analysis_result.nextflow_config or {}
        reasoning = analysis_result.config_reasoning or "基于检测数据的智能分析"
        
        print(f"✅ 配置生成完成")
        
        # 合并配置参数
        final_config = current_config.copy()
        final_config.update(config_params)
        
        return {
            "nextflow_config": final_config,
            "config_reasoning": reasoning,
            "response": f"智能配置分析完成\n\n💡 {reasoning}\n\n🔧 生成了 {len(config_params)} 个配置参数",
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
