from typing import Dict, Any
from ..state import AgentState, PlanResponse
from ..core import get_shared_llm

def create_plan_agent():
    """创建Plan节点的智能规划Agent"""
    llm = get_shared_llm()
    structured_llm = llm.with_structured_output(PlanResponse, method="json_mode")
    return structured_llm

def _build_planning_prompt(state: AgentState, initial_requirements: dict, replan_requirements: dict, is_replanning: bool = False) -> str:
    """构建统一的规划提示词"""
    initial_requirements = initial_requirements or {}
    replan_requirements = replan_requirements or {}
    
    # 构建需求部分
    requirements_section = ""
    if initial_requirements:
        requirements_section += f"\n初始配置需求: {initial_requirements}"
    if replan_requirements:
        requirements_section += f"\n重新规划需求: {replan_requirements}"
    
    # 基础prompt部分
    base_prompt = f"""你是RNA-seq分析{'重新' if is_replanning else ''}规划专家。请基于现有配置状态和用户需求，生成最优化的检测任务计划。

请以JSON格式返回规划结果。

当前nextflow_config: {state.nextflow_config}{requirements_section}

**综合需求处理策略**:
- 如果同时存在初始需求和重新规划需求，优先满足重新规划需求
- 综合考虑所有需求，生成最合适的检测计划
- 确保检测任务能够满足用户的最终配置要求"""

    # 根据是否重新规划添加特定内容
    if is_replanning:
        completed_tasks = list(getattr(state, 'query_results', {}).keys())
        specific_section = f"""
已完成的检测任务: {completed_tasks}

**重新规划策略**:
1. 保留有效的检测结果，避免重复检测
2. 根据综合的需求分析，调整检测策略
3. 优化检测任务的顺序和范围

**综合需求优先处理**:
- 如果任一需求指定新基因组，必须重新执行verify_genome_setup
- 如果任一需求指定新工具，必须重新检测相应工具可用性
- 基于综合需求重新评估必要的检测任务
- 重新规划需求的优先级高于原始需求"""
    else:
        specific_section = f"""

**综合需求处理**:
- 如果任一需求提到特定基因组(如hg38, mm10)，必须执行verify_genome_setup验证可用性
- 如果任一需求提到工具选择，应执行相应的工具可用性检测
- 重新规划需求优先级更高，如有冲突以重新规划需求为准

**智能跳过规则**:
- 如果已配置完整工具链(qc_tool, align_tool, quant_tool) → 可跳过相应工具检测
- verify_genome_setup必须执行，确保基因组文件完整性
- analyze_fastq_data必须执行，这是RNA-seq分析的基础步骤"""

    task_section = """

**可用检测任务**:
1. analyze_fastq_data - FASTQ数据分析和样本配对检测
2. assess_system_readiness - 系统资源和环境准备度评估  
3. verify_genome_setup - 基因组设置和文件完整性验证
4. check_fastp_availability - 检测fastp工具可用性(初次必须)
5. check_star_availability - 检测STAR工具可用性（初次必须）
6. check_featurecounts_availability - 检测featureCounts工具可用性（初次必须）

请返回JSON格式:
- plan: 检测任务列表"""

    return base_prompt + specific_section + task_section

async def plan_node(state: AgentState) -> Dict[str, Any]:
    """增强的Plan节点 - 支持初次规划和重新规划"""
    plan_agent = create_plan_agent()
    
    # 分别获取两种需求
    initial_requirements = getattr(state, 'user_requirements', {})
    replan_requirements = getattr(state, 'replan_requirements', {})
    
    # 判断是否为重新规划
    is_replanning = bool(replan_requirements)
    
    print(f"{'🔄 检测到重新规划请求' if is_replanning else '🎆 初次规划，生成检测计划'}...")
    if initial_requirements:
        print(f"📝 初始配置需求: {initial_requirements}")
    if replan_requirements:
        print(f"🔄 重新规划需求: {replan_requirements}")
    
    # 统一使用一个prompt构建函数
    planning_prompt = _build_planning_prompt(state, initial_requirements, replan_requirements, is_replanning)
    
    try:
        plan_response = await plan_agent.ainvoke([{"role": "user", "content": planning_prompt}])
        detection_plan = plan_response.plan or []
        
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
🔄 **执行策略:** {"基于之前的检测结果智能优化" if is_replanning else "基于现有配置智能优化"}，避免重复检测
{f"📝 **初始配置需求:** {initial_requirements}" if initial_requirements else ""}
{f"🔄 **重新规划需求:** {replan_requirements}" if replan_requirements else ""}

💡 开始执行检测任务..."""
    
    return {
        "plan": detection_plan,
        "response": response_message,
        "status": "plan",
        "user_requirements": initial_requirements,  # 保持初始需求
        "replan_requirements": replan_requirements  # 传递重新规划需求
    }