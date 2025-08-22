from typing import Dict, Any
from ..state import AgentState, PlanResponse
from ..core import get_shared_llm

def create_plan_agent():
    """创建Plan节点的智能规划Agent"""
    llm = get_shared_llm()
    structured_llm = llm.with_structured_output(PlanResponse, method="json_mode")
    return structured_llm

async def plan_node(state: AgentState) -> Dict[str, Any]:
    """增强的Plan节点 - 支持初次规划和重新规划"""
    plan_agent = create_plan_agent()
    
    # 检查用户是否有具体需求
    user_requirements = ""
    if state.messages and len(state.messages) > 0:
        latest_message = state.messages[-1]
        if isinstance(latest_message, dict):
            content = latest_message.get("content", "")
        else:
            content = getattr(latest_message, "content", "")
        
        # 提取用户的具体需求（去掉命令部分）
        if "/replan" in content.lower():
            user_requirements = content.lower().replace("/replan", "").strip()
        elif "/plan" in content.lower():
            user_requirements = content.lower().replace("/plan", "").strip()
        else:
            user_requirements = content.strip()
    
    # 检查是否为重新规划
    is_replanning = bool(state.query_results) or "/replan" in (state.messages[-1].get("content", "") if state.messages else "")
    
    if is_replanning:
        print("🔄 检测到重新规划请求，整合用户需求...")
        if user_requirements:
            print(f"📝 用户需求: {user_requirements}")
        planning_prompt = _build_replanning_prompt(state, user_requirements)
    else:
        print("🎆 初次规划，生成检测计划...")
        if user_requirements:
            print(f"📝 用户需求: {user_requirements}")
        planning_prompt = _build_initial_planning_prompt(state, user_requirements)
    
    try:
        plan_response = plan_agent.invoke(planning_prompt)
        detection_plan = plan_response.plan or []
        analysis_intent = plan_response.analysis_intent or "RNA-seq标准分析"
        
        if not detection_plan:
            raise Exception("LLM未生成有效的检测计划")
            
    except Exception as e:
        print(f"❌ LLM规划失败: {e}")
        return {
            "plan": [],
            "response": f"❌ 规划失败: {str(e)}\n\n💡 请重新尝试或检查系统配置",
            "status": "normal"  # 路由回normal模式
        }
    
    response_message = f"""🎆 **{"重新" if is_replanning else ""}智能分析计划制定完成**

📋 **优化检测计划:** {len(detection_plan)} 个任务
🎭 **分析目标:** {analysis_intent}
{f"📝 **用户需求:** {user_requirements}" if user_requirements else ""}
🔄 **执行策略:** {"基于之前的检测结果智能优化" if is_replanning else "基于现有配置智能优化"}，避免重复检测

💡 开始执行检测任务..."""
    
    return {
        "plan": detection_plan,
        "response": response_message,
        "status": "plan"
    }

def _build_initial_planning_prompt(state: AgentState, user_requirements: str = "") -> str:
    """构建初次规划提示词"""
    user_section = f"\n用户特殊需求: {user_requirements}" if user_requirements else ""
    
    return f"""你是RNA-seq分析规划专家。请基于现有配置状态，生成最优化的检测任务计划。

当前nextflow_config: {state.nextflow_config}{user_section}

**可用检测任务**:
1. analyze_fastq_data - FASTQ数据分析和样本配对检测
2. assess_system_readiness - 系统资源和环境准备度评估  
3. verify_genome_setup - 基因组设置和文件完整性验证
4. check_fastp_availability - 检测fastp工具可用性
5. check_star_availability - 检测STAR工具可用性
6. check_featurecounts_availability - 检测featureCounts工具可用性

**智能跳过规则**:
- 如果已配置genome_version且species → 可跳过verify_genome_setup
- 如果已配置完整工具链(qc_tool, align_tool, quant_tool) → 可跳过assess_system_readiness

**用户需求处理**:
- 如果用户提到特定基因组(如hg38, mm10)，应重点执行verify_genome_setup
- 如果用户提到工具选择，应执行相应的工具可用性检测

请返回JSON格式:
{{
  "plan": [检测任务列表],
  "analysis_intent": "分析目标描述"
}}"""

def _build_replanning_prompt(state: AgentState, user_requirements: str = "") -> str:
    """构建重新规划提示词"""
    completed_tasks = list(state.query_results.keys()) if state.query_results else []
    user_section = f"\n用户新需求: {user_requirements}" if user_requirements else ""
    
    return f"""你是RNA-seq分析重新规划专家。用户要求重新规划，请基于已有的检测结果优化计划。

当前nextflow_config: {state.nextflow_config}
已完成的检测任务: {completed_tasks}{user_section}

**重新规划策略**:
1. 保留有效的检测结果，避免重复检测
2. 根据用户的新需求，调整检测策略
3. 优化检测任务的顺序和范围

**用户需求优先处理**:
- 如果用户指定新基因组，必须重新执行verify_genome_setup
- 如果用户指定新工具，必须重新检测相应工具可用性
- 基于用户需求重新评估必要的检测任务

**可选的重新检测任务**:
- analyze_fastq_data (如果数据源发生变化)
- assess_system_readiness (如果需要重新评估资源)
- verify_genome_setup (如果基因组配置有问题)
- check_fastp_availability (如果需要重新检测fastp)
- check_star_availability (如果需要重新检测STAR)
- check_featurecounts_availability (如果需要重新检测featureCounts)

请返回JSON格式:
{{
  "plan": [优化后的检测任务列表],
  "analysis_intent": "重新规划的目标和理由"
}}"""