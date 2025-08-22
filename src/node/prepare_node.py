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
    """准备节点 - 基于检测数据让LLM直接生成配置参数"""
    print(f"⚙️ 开始智能配置分析...")
    
    # 直接从state获取检测数据
    detection_results = state.query_results or {}
    current_config = state.nextflow_config or {}
    
    if not detection_results:
        print("⚠️ 未检测到任何数据，无法生成配置")
        return {
            "nextflow_config": current_config,
            "config_reasoning": "未获取到检测数据，保持现有配置",
            "response": "⚠️ 缺少检测数据，无法进行智能配置分析",
            "status": "error"
        }
    
    # 使用LLM直接生成配置
    prepare_agent = create_prepare_agent()
    
    try:
        print("🧠 LLM正在分析检测数据并生成配置...")
        
        # 构建简单的系统消息
        system_message = """你是RNA-seq分析配置专家。基于检测数据直接生成Nextflow配置参数和分析理由。

重要：请返回有效的JSON格式，包含以下字段：
- nextflow_config: 配置参数对象
- config_reasoning: 配置理由字符串

示例：
{
  "nextflow_config": {"genome_version": "hg38", "qc_tool": "fastp", "threads": 8},
  "config_reasoning": "基于检测数据选择标准工具链"
}"""
        
        # 直接使用原始检测数据，不需要格式化
        user_message = f"""请基于以下检测数据生成配置参数：

检测数据 (JSON格式):
{json.dumps(detection_results, indent=2, ensure_ascii=False)}

当前配置: {current_config}

请分析检测数据并生成合适的Nextflow配置参数。"""
        
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
