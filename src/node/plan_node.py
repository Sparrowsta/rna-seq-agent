from typing import Dict, Any
from ..state import AgentState, PlanResponse
from ..core import get_shared_llm
from ..prompts import PLAN_NODE_PROMPT
from langgraph.prebuilt import create_react_agent

def create_plan_agent(requirements_context: str = ""):
    """创建Plan节点的智能规划Agent
    
    Args:
        requirements_context: 用户需求和状态的上下文信息
    """
    llm = get_shared_llm()
    
    # 构建完整的prompt，结合基础prompt和上下文
    if requirements_context:
        full_prompt = PLAN_NODE_PROMPT + "\n\n## 当前分析上下文\n" + requirements_context
    else:
        full_prompt = PLAN_NODE_PROMPT
    
    # 使用create_react_agent但不提供tools，纯推理模式
    agent = create_react_agent(
        model=llm,
        tools=[],  # 空工具列表，纯推理
        prompt=full_prompt,  # 使用动态构建的完整prompt
        response_format=PlanResponse
    )
    return agent

async def plan_node(state: AgentState) -> Dict[str, Any]:
    """增强的Plan节点 - 支持并行任务组规划"""
    
    # 分别获取两种需求
    initial_requirements = getattr(state, 'user_requirements', {})
    replan_requirements = getattr(state, 'replan_requirements', {})
    
    # 判断是否为重新规划
    is_replanning = bool(replan_requirements)
    
    # 构建上下文信息
    context_parts = []
    if initial_requirements:
        context_parts.append(f"初始配置需求: {initial_requirements}")
    if replan_requirements:
        context_parts.append(f"重新规划需求: {replan_requirements}")
    
    # 添加规划类型说明
    if is_replanning:
        context_parts.append("这是一次重新规划，请基于新的需求调整检测计划")
    else:
        context_parts.append("这是初始规划，请生成完整的检测任务计划")
        
    # 添加当前配置状态
    if hasattr(state, 'nextflow_config') and state.nextflow_config:
        context_parts.append(f"当前配置状态: {state.nextflow_config}")
    
    requirements_context = "\n".join(context_parts) if context_parts else ""
    
    try:
        # 将上下文传递给create_plan_agent
        agent_executor = create_plan_agent(requirements_context)
        
        # 简单的用户消息
        user_message = "请生成并行检测任务计划"
        messages_input = {"messages": [{"role": "user", "content": user_message}]}
        
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")
        
        if structured_response:
            # 提取并行任务组信息
            task_groups = structured_response.plan or []
            group_descriptions = structured_response.group_descriptions or []
            execution_strategy = structured_response.execution_strategy or "parallel"
        
        if not task_groups:
            raise Exception("未生成有效的检测计划")
            
    except Exception as e:
        print(f"❌ LLM规划失败: {e}")
        return {
            "plan": [],
            "group_descriptions": [],
            "execution_strategy": "sequential",
            "error": f"规划失败: {str(e)}"
        }
    
    # 平铺任务组为检测任务列表
    detection_tasks = []
    for group in task_groups:
        detection_tasks.extend(group)
    
    print(f"📋 生成检测计划: {len(task_groups)}组并行任务，共{len(detection_tasks)}个检测项")
    for i, (group, desc) in enumerate(zip(task_groups, group_descriptions)):
        print(f"  组{i+1}: {desc} -> {group}")
    
    return {
        "plan": task_groups,
        "group_descriptions": group_descriptions, 
        "execution_strategy": execution_strategy,
        "detection_tasks": detection_tasks,
        "is_replanning": is_replanning
    }