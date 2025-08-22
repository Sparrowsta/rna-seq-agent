from typing import Dict, Any
from ..state import AgentState, ModifyResponse
from ..core import get_shared_llm

def create_modify_agent():
    """创建Modify节点的配置修改Agent"""
    llm = get_shared_llm()
    structured_llm = llm.with_structured_output(ModifyResponse, method="json_mode")
    return structured_llm

async def modify_node(state: AgentState) -> Dict[str, Any]:
    """Modify节点 - 接受用户修改请求并直接返回修改后的nextflow配置"""
    # 从messages中获取最新的用户输入
    user_input = ""
    if state.messages:
        # 获取最后一条用户消息
        for msg in reversed(state.messages):
            if isinstance(msg, dict) and msg.get('role') == 'user':
                user_input = msg.get('content', '')
                break
            elif isinstance(msg, str):
                user_input = msg
                break
    
    current_config = state.nextflow_config or {}
    
    modify_agent = create_modify_agent()
    
    modify_prompt = f"""你是RNA-seq配置修改专家。用户想要修改当前的分析配置，请根据用户请求修改nextflow参数。

当前配置: {current_config}
用户修改请求: {user_input}

**可修改的参数:**
- genome_version: 基因组版本 (如: hg38, hg19, mm10, mm39等)
- species: 物种 (如: human, mouse等)
- qc_tool: 质控工具 (fastp, cutadapt)
- align_tool: 比对工具 (star, hisat2)  
- quant_tool: 定量工具 (featurecounts, htseq)

**修改原则:**
- 只修改用户明确提到的参数
- 保持其他参数不变
- 确保参数组合的兼容性
- 如果用户请求不清晰，保持当前配置不变

请返回JSON格式，包含:
- modified_config: 修改后的完整配置
- modification_summary: 修改内容的简要说明

基于用户请求修改配置参数。"""
    
    try:
        modify_response = modify_agent.invoke(modify_prompt)
        modified_config = modify_response.modified_config or current_config
        modification_summary = modify_response.modification_summary or "配置未发生变化"
    except Exception:
        # 如果LLM调用失败，保持原配置不变
        modified_config = current_config
        modification_summary = "配置修改失败，保持原配置"
    
    response_message = f"""🛠️ **配置修改完成**

📝 **修改说明:** {modification_summary}
⚙️ **新配置已生成，等待确认**

请确认修改后的配置是否正确。"""
    
    return {
        "nextflow_config": modified_config,
        "config_reasoning": modification_summary,
        "response": response_message,
        "status": "confirm"
    }