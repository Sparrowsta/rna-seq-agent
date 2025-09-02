import json
from typing import Dict, Any
from ..state import AgentState, PrepareResponse
from ..core import get_shared_llm
from ..prompts import PREPARE_NODE_PROMPT
from langgraph.prebuilt import create_react_agent

def create_prepare_agent(detection_context: str = ""):
    """创建Prepare节点的智能配置Agent
    
    Args:
        detection_context: 检测数据和用户需求的上下文信息
    """
    llm = get_shared_llm()
    
    # 构建完整的prompt，结合基础prompt和上下文
    if detection_context:
        full_prompt = PREPARE_NODE_PROMPT + "\n\n## 当前分析数据\n" + detection_context
    else:
        full_prompt = PREPARE_NODE_PROMPT
    
    # 使用create_react_agent但不提供tools，纯推理模式
    agent = create_react_agent(
        model=llm,
        tools=[],  # 空工具列表，纯推理
        prompt=full_prompt,  # 使用动态构建的完整prompt
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
        
        # 构建需求上下文
        context_parts = []
        if initial_requirements:
            context_parts.append(f"初始配置需求: {initial_requirements}")
        if replan_requirements:
            context_parts.append(f"重新规划需求: {replan_requirements} (优先级更高)")
        
        # 添加检测数据
        context_parts.append(f"=== 📊 系统检测数据 ===")
        context_parts.append(json.dumps(detection_results, indent=2, ensure_ascii=False))
        
        # 添加当前配置状态
        context_parts.append(f"=== ⚙️ 当前配置状态 ===")
        context_parts.append(json.dumps(current_config, indent=2, ensure_ascii=False))
        
        detection_context = "\n".join(context_parts)
        
        # 将上下文传递给create_prepare_agent
        agent_executor = create_prepare_agent(detection_context)
        
        # 简单的用户消息
        user_message = "请基于用户需求和检测数据生成最优配置"
        messages_input = {"messages": [{"role": "user", "content": user_message}]}
        
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")

        # 检查LLM响应并提取结果（严格遵循 PrepareResponse 的字段：nextflow_config/resource_config/config_reasoning）
        if structured_response:
            reasoning = structured_response.config_reasoning or "基于用户需求和检测数据的智能分析"

            nextflow_cfg = structured_response.nextflow_config or {}
            resource_params = structured_response.resource_config or {}

            print(f"✅ 配置生成完成，严格遵循用户需求")

            # 合并配置参数（新配置优先）
            final_config = current_config.copy()
            final_config.update(nextflow_cfg)

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
                "resource_config": resource_params,  # 显式传递资源配置，供 execute_node 生成 nextflow.config
                "config_reasoning": reasoning,
                "response": f"智能配置分析完成{user_satisfaction_note}\n\n💡 {reasoning}",
                "status": "confirm"
            }
        else:
            raise Exception("Agent未返回预期的结构化响应")
        
    except Exception as e:
        print(f"❌ LLM分析失败: {str(e)}")
        return {
            "nextflow_config": current_config,
            "resource_config": {},
            "config_reasoning": f"LLM分析失败: {str(e)}",
            "response": f"❌ 配置生成失败: {str(e)}",
            "status": "error"
        }
