from typing import Dict, Any
from ..state import AgentState, PlanResponse
from ..core import get_shared_llm
from langgraph.prebuilt import create_react_agent

def create_plan_agent():
    """创建Plan节点的智能规划Agent"""
    llm = get_shared_llm()
    
    # 使用create_react_agent但不提供tools，纯推理模式
    agent = create_react_agent(
        model=llm,
        tools=[],  # 空工具列表，纯推理
        prompt="你是RNA-seq分析规划专家。请基于用户提供的详细需求和配置状态，进行智能推理并生成最优化的检测任务计划。",
        response_format=PlanResponse
    )
    return agent

def _build_planning_prompt(state: AgentState, initial_requirements: dict, replan_requirements: dict, is_replanning: bool = False) -> str:
    """构建并行任务规划的统一提示词"""
    initial_requirements = initial_requirements or {}
    replan_requirements = replan_requirements or {}
    
    # 构建需求部分
    requirements_section = ""
    if initial_requirements:
        requirements_section += f"\n初始配置需求: {initial_requirements}"
    if replan_requirements:
        requirements_section += f"\n重新规划需求: {replan_requirements}"
    
    # 基础prompt部分
    base_prompt = f"""你是RNA-seq分析{'重新' if is_replanning else ''}规划专家。请基于现有配置状态和用户需求，生成最优化的并行检测任务计划。

请以JSON格式返回规划结果。

当前nextflow_config: {state.nextflow_config}{requirements_section}

**并行执行策略**:
- 将所有检测任务分为6个独立的串行组
- 每个串行组可以并行执行，组内任务按顺序执行
- 这样设计便于将来在各组内添加有依赖关系的任务"""

    # 根据是否重新规划添加特定内容
    if is_replanning:
        completed_tasks = list(getattr(state, 'query_results', {}).keys())
        specific_section = f"""
已完成的检测任务: {completed_tasks}

**重新规划策略**:
1. 保留有效的检测结果，避免重复检测
2. 根据综合的需求分析，调整检测策略
3. 优化检测任务的顺序和分组

**综合需求优先处理**:
- 如果任一需求指定新基因组，必须重新执行verify_genome_setup
- 如果任一需求指定新工具，必须重新检测相应工具可用性
- 基于综合需求重新评估必要的检测任务
- 重新规划需求的优先级高于原始需求"""
    else:
        specific_section = f"""

**并行分组原则**:
- 每个检测任务都是独立的，可以同时执行
- 分为6个独立串行组，最大化并行度
- 组内预留扩展空间，便于将来添加依赖任务

**智能跳过规则**:
- 如果已配置完整工具链(qc_tool, align_tool, quant_tool) → 可跳过相应工具检测
- verify_genome_setup必须执行，确保基因组文件完整性
- analyze_fastq_data必须执行，这是RNA-seq分析的基础步骤"""

    task_section = """

**7个独立并行组设计**:
1. 数据检测组: ["analyze_fastq_data"] - FASTQ数据分析和样本配对检测
2. 系统检测组: ["assess_system_readiness"] - 系统资源和环境准备度评估  
3. QC工具组: ["check_fastp_availability"] - 检测fastp工具可用性
4. STAR比对组: ["check_star_availability"] - 检测STAR工具可用性
5. HISAT2比对组: ["check_hisat2_availability"] - 检测HISAT2工具可用性
6. 定量工具组: ["check_featurecounts_availability"] - 检测featureCounts工具可用性
7. 基因组配置组: ["verify_genome_setup"] - 基因组设置和文件完整性验证

请返回JSON格式:
- plan: [[组1任务列表], [组2任务列表], ...]
- group_descriptions: ["组1描述", "组2描述", ...]
- execution_strategy: "parallel"""

    return base_prompt + specific_section + task_section

async def plan_node(state: AgentState) -> Dict[str, Any]:
    """增强的Plan节点 - 支持并行任务组规划"""
    
    # 分别获取两种需求
    initial_requirements = getattr(state, 'user_requirements', {})
    replan_requirements = getattr(state, 'replan_requirements', {})
    
    # 判断是否为重新规划
    is_replanning = bool(replan_requirements)
    
    # 统一使用一个prompt构建函数
    planning_prompt = _build_planning_prompt(state, initial_requirements, replan_requirements, is_replanning)
    
    try:
        agent_executor = create_plan_agent()
        messages_input = {"messages": [{"role": "user", "content": planning_prompt}]}
        
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")
        
        if structured_response:
            # 提取并行任务组信息
            task_groups = structured_response.plan or []
            group_descriptions = structured_response.group_descriptions or []
            execution_strategy = structured_response.execution_strategy or "parallel"
                
            # 如果没有分组信息，使用默认7组分组
            if not task_groups:
                task_groups = [
                    ["analyze_fastq_data"],
                    ["assess_system_readiness"],
                    ["check_fastp_availability"],
                    ["check_star_availability"],
                    ["check_hisat2_availability"],
                    ["check_featurecounts_availability"],
                    ["verify_genome_setup"]
                ]
                group_descriptions = [
                    "数据检测组",
                    "系统检测组", 
                    "QC工具组",
                    "STAR比对组",
                    "HISAT2比对组",
                    "定量工具组",
                    "基因组配置组"
                ]
        else:
            # 如果没有结构化响应，使用默认计划
            task_groups = [
                ["analyze_fastq_data"],
                ["assess_system_readiness"],
                ["check_fastp_availability"],
                ["check_star_availability"],
                ["check_hisat2_availability"],
                ["check_featurecounts_availability"],
                ["verify_genome_setup"]
            ]
            group_descriptions = [
                "数据检测组", 
                "系统检测组", 
                "QC工具组",
                "STAR比对组",
                "HISAT2比对组", 
                "定量工具组",
                "基因组配置组"
            ]
            execution_strategy = "parallel"
        
        if not task_groups:
            raise Exception("未生成有效的检测计划")
            
    except Exception as e:
        print(f"❌ LLM规划失败: {e}")
        return {
            "plan": [],
            "group_descriptions": [],
            "execution_strategy": "sequential",
            "response": f"❌ 规划失败: {str(e)}\n\n💡 请重新尝试或检查系统配置",
            "status": "normal"  # 路由回normal模式
        }
    
    response_message = f"""🎆 **{"重新" if is_replanning else ""}智能并行分析计划制定完成**

📋 **并行任务组:** {len(task_groups)} 个独立组，总计 {sum(len(group) for group in task_groups)} 个任务
🚀 **执行策略:** {execution_strategy} - {"所有组同时执行，最大化效率" if execution_strategy == "parallel" else "顺序执行"}
🔄 **优化策略:** {"基于之前的检测结果智能优化" if is_replanning else "基于现有配置智能优化"}
{f"📝 **初始配置需求:** {initial_requirements}" if initial_requirements else ""}
{f"🔄 **重新规划需求:** {replan_requirements}" if replan_requirements else ""}

🎯 **任务分组详情:**
{chr(10).join([f"  组{i+1}: {desc} -> {group}" for i, (desc, group) in enumerate(zip(group_descriptions, task_groups))])}

💡 开始执行并行检测任务..."""
    
    return {
        "plan": task_groups,  # 直接存储任务组
        "group_descriptions": group_descriptions,
        "execution_strategy": execution_strategy,
        "response": response_message,
        "status": "detect",
        "user_requirements": initial_requirements,  # 保持初始需求
        "replan_requirements": replan_requirements  # 传递重新规划需求
    }