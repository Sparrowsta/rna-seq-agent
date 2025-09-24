"""
FeatureCounts节点 - 用于执行基因定量分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FeaturecountsResponse
from ..core import get_shared_llm
from ..prompts import FEATURECOUNTS_OPTIMIZATION_PROMPT
from ..tools import (
    run_nextflow_featurecounts,
    parse_featurecounts_metrics,
    scan_genome_files,
    extract_genome_paths
)
from ..route_decider import decide_next_action_featurecounts
from ..logging_bootstrap import get_logger, log_llm_preview
import json
from datetime import datetime
from pathlib import Path


def create_featurecounts_agent():
    """创建FeatureCounts节点的React Agent"""
    llm = get_shared_llm()
    
    # 绑定FeatureCounts相关工具
    tools = [
        scan_genome_files,
        run_nextflow_featurecounts, 
        parse_featurecounts_metrics
    ]
    
    # 创建React Agent，使用专业的FeatureCounts优化Prompt
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=FEATURECOUNTS_OPTIMIZATION_PROMPT,
        response_format=FeaturecountsResponse
    )
    return agent

def append_featurecounts_optimization_history(state: AgentState, optimization_params: Dict[str, Any], 
                                           suggestions: str, results: Dict[str, Any]) -> None:
    """追加FeatureCounts优化历史记录，保持最近5次记录"""
    history_entry = {
        "timestamp": datetime.now().isoformat(),
        "execution_id": f"featurecounts_run_{len(state.featurecounts_optimization_history) + 1}",
        "optimization_params": optimization_params or {},
        "optimization_suggestions": suggestions or "",
        "execution_results": results or {}
    }
    
    # 追加新记录
    state.featurecounts_optimization_history.append(history_entry)
    
    # 保持最近5次记录
    if len(state.featurecounts_optimization_history) > 5:
        state.featurecounts_optimization_history = state.featurecounts_optimization_history[-5:]
    
    logger = get_logger("rna.nodes.featurecounts")
    logger.info(f"[FEATURECOUNTS] 已追加优化历史记录，当前保存{len(state.featurecounts_optimization_history)}次历史")


async def featurecounts_node(state: AgentState) -> Dict[str, Any]:
    """
    FeatureCounts节点实现 - 执行基因定量并生成优化建议
    
    功能：
    - 执行基因定量
    - 生成表达矩阵  
    - 更新状态信息
    - 生成optimization_params供路由决策器使用
    """
    logger = get_logger("rna.nodes.featurecounts")
    logger.info("FeatureCounts定量节点开始执行")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "featurecounts" not in completed_steps:
        completed_steps.append("featurecounts")
    
    # 允许基于 STAR 或 HISAT2 的比对结果进行定量
    has_star = bool(getattr(state, 'star_results', {}) or {}) and bool(state.star_results.get("success"))
    has_hisat2 = bool(getattr(state, 'hisat2_results', {}) or {}) and bool(state.hisat2_results.get("success"))
    if not (has_star or has_hisat2):
        # 依赖失败时设置返回上下文
        state.return_source = "featurecounts"
        state.return_reason = "failed"

        return {
            "success": False,
            "status": "featurecounts_failed",
            "response": "❌ FeatureCounts执行失败：缺少有效的比对结果（STAR/HISAT2），请先完成比对",
            "current_step": "featurecounts",
            "completed_steps": completed_steps,
            "featurecounts_results": {
                "success": False,
                "status": "failed",
                "error": "比对结果不可用或未成功"
            }
        }
    
    try:
        logger.info("[FEATURECOUNTS] 使用FeatureCounts Agent进行定量分析和优化...")
        
        # 统一调用FeatureCounts Agent执行定量和优化分析
        fc_response = await _call_featurecounts_optimization_agent(state)
        
        # 立即更新执行参数和优化建议
        optimized_params = fc_response.featurecounts_params
        optimization_reasoning = fc_response.featurecounts_optimization_suggestions
        optimization_params_changes = fc_response.featurecounts_optimization_params
        
        # 处理执行结果（以工具真实返回为准，默认失败避免"空成功"）
        try:
            agent_results = dict(getattr(fc_response, 'featurecounts_results', {}) or {})
        except Exception:
            agent_results = {}
        success_flag = bool(agent_results.get("success", False))
        status_text = agent_results.get("status", "success" if success_flag else "failed")
        fc_results = {**agent_results, "success": success_flag, "status": status_text}

        # 生成响应信息
        optimization_count = len(optimization_params_changes or {})
        if success_flag:
            if optimization_count > 0:
                response = (
                    f"✅ FeatureCounts定量完成\n- 定量状态: 成功完成\n- 优化分析: 生成了{optimization_count}个优化建议\n\n"
                    f"⚡ 优化详情: {optimization_reasoning}"
                )
            else:
                response = (
                    "✅ FeatureCounts定量完成\n\n"
                    "🚀 执行详情: 已完成基因定量，当前参数配置已是最优"
                )
        else:
            error_msg = fc_results.get("error") or fc_results.get("message") or "FeatureCounts执行未产生有效输出"
            response = f"❌ FeatureCounts执行失败：{error_msg}"

        logger.info(f"[FEATURECOUNTS] FeatureCounts执行完成，生成{optimization_count}个优化参数")

        # 更新状态以便路由决策读取最新结果
        state.featurecounts_results = fc_results

        # 根据路由决策器结果设置返回上下文
        next_action = decide_next_action_featurecounts(state)
        if next_action == "return_confirm":
            state.return_source = "featurecounts"
            if not success_flag:
                state.return_reason = "failed"
            elif state.execution_mode == 'batch_optimize':
                state.return_reason = "batch_collect"  # Batch模式特殊处理
            else:
                state.return_reason = "step_confirm"

        # 构建成功结果
        result = {
            "success": success_flag,
            "status": "featurecounts_completed",
            "current_step": "featurecounts",
            "completed_steps": completed_steps,
            "response": response,
            "featurecounts_params": optimized_params,
            "featurecounts_optimization_suggestions": optimization_reasoning,
            "featurecounts_optimization_params": optimization_params_changes,
            "featurecounts_results": fc_results,
        }

        # 追加优化历史记录
        append_featurecounts_optimization_history(
            state=state,
            optimization_params=optimization_params_changes,
            suggestions=optimization_reasoning,
            results=fc_results
        )

        logger.info(f"[FEATURECOUNTS] FeatureCounts执行完成，生成{optimization_count}个优化参数")
        return result
            
    except Exception as e:
        logger.error(f"FeatureCounts节点执行失败: {str(e)}", exc_info=True)

        # 失败时设置返回上下文
        state.return_source = "featurecounts"
        state.return_reason = "failed"

        return {
            "success": False,
            "status": "featurecounts_failed",
            "response": f"❌ FeatureCounts定量执行失败: {str(e)}",
            "current_step": "featurecounts",
            "completed_steps": completed_steps,
            "featurecounts_results": {
                "success": False,
                "status": "failed",
                "error": str(e)
            },
        }


async def _call_featurecounts_optimization_agent(state: AgentState) -> FeaturecountsResponse:
    """调用FeatureCounts优化Agent进行智能参数优化"""
    
    logger = get_logger("rna.nodes.featurecounts")
    
    # 获取基因组配置信息（从detect节点的query_results中获取）
    # 动态提取基因组路径信息
    genome_paths = extract_genome_paths(state)
    genome_version = state.nextflow_config.get("genome_version")
    
    # 提取GTF路径用于工具调用
    gtf_path = genome_paths.get("gtf_path", "")
    
    # 提取资源配置片段（仅 FeatureCounts）
    featurecounts_resource_config = state.resource_config.get("featurecounts") if state.resource_config else {}

    user_context = {
        "execution_mode": state.execution_mode,
        "genome_config": {
            "genome_version": genome_version,
            "paired_end": state.nextflow_config.get("paired_end")
        },
        "gtf_path": gtf_path,  # 简化传递，只提供GTF路径
        "featurecounts_resource_config": featurecounts_resource_config,
        "current_featurecounts_params": state.featurecounts_params,
        "star_results": getattr(state, 'star_results', {}),
        "hisat2_results": getattr(state, 'hisat2_results', {}),
        "optimization_history": state.featurecounts_optimization_history
    }
    
    user_prompt = json.dumps(user_context, ensure_ascii=False, indent=2)
    
    # 创建并调用Agent
    agent_executor = create_featurecounts_agent()
    
    # 构建消息格式
    messages = [
        {"role": "user", "content": user_prompt}
    ]
    
    result = await agent_executor.ainvoke({"messages": messages})
    
    # 提取结构化响应
    structured_response = result.get("structured_response")
    try:
        if structured_response:
            log_llm_preview(logger, "featurecounts", structured_response)
        else:
            log_llm_preview(logger, "featurecounts.raw", {"keys": list(result.keys())[:10]})
    except Exception:
        pass

    # 定义最小校验：必须含有 results_dir 和关键产物文件存在
    def _is_valid_featurecounts_results(res: FeaturecountsResponse) -> bool:
        try:
            fcr = getattr(res, 'featurecounts_results', {}) or {}
            if not fcr.get('success'):
                return False
            results_dir = fcr.get('results_dir') or fcr.get('results_directory')
            if not results_dir:
                return False

            # 检查关键输出文件是否存在
            output = fcr.get('output', {})
            required_files = [
                output.get('counts_file'),
                output.get('summary_file')
            ]

            # 检查每样本输出文件（如果存在）
            per_sample_outputs = output.get('per_sample_outputs', [])
            for item in per_sample_outputs:
                if item.get('counts_file'):
                    required_files.append(item.get('counts_file'))

            missing_paths = []
            for file_path in required_files:
                if file_path and not Path(file_path).exists():
                    missing_paths.append(file_path)

            return len(missing_paths) == 0
        except Exception:
            return False

    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")

    if not _is_valid_featurecounts_results(structured_response):
        raise ValueError("Agent返回的结果无效或缺少必要产物")

    return structured_response
