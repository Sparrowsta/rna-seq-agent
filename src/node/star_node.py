"""
STAR节点 - 用于执行STAR比对分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, StarResponse
from ..core import get_shared_llm
from ..prompts import STAR_OPTIMIZATION_PROMPT
from ..tools import (
    download_genome_assets,
    build_star_index,
    run_nextflow_star,
    parse_star_metrics,
    extract_genome_paths
)
from ..route_decider import decide_next_action_star
from ..logging_bootstrap import get_logger, log_llm_preview
from pathlib import Path
import json
from datetime import datetime

logger = get_logger("rna.nodes.star")


def create_star_agent():
    """创建STAR节点的React Agent"""
    llm = get_shared_llm()

    # 绑定STAR相关工具
    tools = [
        download_genome_assets,
        build_star_index,
        run_nextflow_star,
        parse_star_metrics,
    ]

    # 创建使用工具的React Agent
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=STAR_OPTIMIZATION_PROMPT,
        response_format=StarResponse,
    )
    return agent


def append_star_optimization_history(state: AgentState, optimization_params: Dict[str, Any],
                                   suggestions: str, results: Dict[str, Any]) -> None:
    """追加STAR优化历史记录，保持最近5次记录"""
    history_entry = {
        "timestamp": datetime.now().isoformat(),
        "execution_id": f"star_run_{len(state.star_optimization_history) + 1}",
        "optimization_params": optimization_params or {},
        "optimization_suggestions": suggestions or "",
        "execution_results": results or {}
    }

    # 追加新记录
    state.star_optimization_history.append(history_entry)

    # 保持最近5次记录
    if len(state.star_optimization_history) > 5:
        state.star_optimization_history = state.star_optimization_history[-5:]

    logger.info(f"[STAR] 已追加优化历史记录，当前保存{len(state.star_optimization_history)}次历史")
    logger.info(f"[STAR] 最新优化历史记录：{state.star_optimization_history[-1]}")


async def star_node(state: AgentState) -> Dict[str, Any]:
    """
    STAR节点实现 - 执行序列比对并生成优化建议
    
    功能：
    - 执行STAR比对
    - 基于Agent智能优化参数
    - 更新状态信息
    - 生成optimization_params供路由决策器使用
    """
    logger.info("STAR比对节点开始执行...")

    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "star" not in completed_steps:
        completed_steps.append("star")

    # 依赖检查：需要FastP成功结果
    if not state.fastp_results or not state.fastp_results.get("success"):
        # 依赖失败时设置返回上下文
        state.return_source = "star"
        state.return_reason = "failed"

        return {
            "success": False,
            "status": "star_failed",
            "response": "❌ STAR比对失败：缺少有效的FastP质控结果，请先完成FastP质量控制",
            "current_step": "star",
            "completed_steps": completed_steps,
            "star_results": {
                "success": False,
                "status": "failed",
                "error": "missing_fastp_results"
            }
        }

    try:
        # 统一通过Agent执行STAR比对和优化分析
        logger.info("[STAR] 调用Agent执行STAR比对和优化分析...")
        agent_response = await _call_star_optimization_agent(state)

        # 更新执行参数和优化建议
        optimized_params = agent_response.star_params
        optimization_reasoning = agent_response.star_optimization_suggestions
        optimization_params_changes = agent_response.star_optimization_params

        # 处理执行结果（以工具真实返回为准，默认失败避免“空成功”）
        try:
            agent_results = dict(getattr(agent_response, 'star_results', {}) or {})
        except Exception:
            agent_results = {}
        success_flag = bool(agent_results.get("success", False))
        status_text = agent_results.get("status", "success" if success_flag else "failed")
        star_results = {**agent_results, "success": success_flag, "status": status_text}

        # 生成响应信息
        optimization_count = len(optimization_params_changes or {})
        if success_flag:
            if optimization_count > 0:
                response = (
                    f"✅ STAR比对完成\n- 比对状态: 成功完成\n- 优化分析: 生成了{optimization_count}个优化建议\n\n"
                    f"⚡ 优化详情: {optimization_reasoning}"
                )
            else:
                response = (
                    "✅ STAR比对完成\n\n"
                    "🚀 执行详情: 已完成序列比对"
                )
        else:
            error_msg = star_results.get("error") or star_results.get("message") or "STAR执行未产生有效输出"
            response = f"❌ STAR执行失败：{error_msg}"
            
        logger.info(f"[STAR] STAR执行完成，生成{optimization_count}个优化参数")

        # 追加优化历史记录
        append_star_optimization_history(
            state=state,
            optimization_params=optimization_params_changes,
            suggestions=optimization_reasoning,
            results=star_results
        )

        # 更新状态以便路由决策读取最新结果
        state.star_results = star_results

        # 根据路由决策器结果设置返回上下文
        next_action = decide_next_action_star(state)
        if next_action == "return_confirm":
            state.return_source = "star"
            if not success_flag:
                state.return_reason = "failed"
            elif state.execution_mode == 'batch_optimize' and optimization_count > 0:
                state.return_reason = "batch_collect"
            else:
                state.return_reason = "step_confirm"

        # 构建成功结果
        result = {
            "success": success_flag,
            "status": "star_completed",
            "current_step": "star",
            "completed_steps": completed_steps,
            "response": response,
            "star_params": optimized_params,
            "star_optimization_suggestions": optimization_reasoning,
            "star_optimization_params": optimization_params_changes,
            "star_results": star_results,
        }

        return result

    except Exception as e:
        logger.error(f"[STAR] STAR执行失败: {str(e)}")

        # 失败时设置返回上下文
        state.return_source = "star"
        state.return_reason = "failed"

        return {
            "success": False,
            "status": "star_failed",
            "response": f"❌ STAR执行失败: {str(e)}",
            "current_step": "star",
            "completed_steps": completed_steps,
            "star_results": {
                "success": False,
                "status": "failed",
                "error": str(e)
            }
        }


async def _call_star_optimization_agent(state: AgentState) -> StarResponse:
    """调用STAR优化Agent进行智能参数优化"""

    # 动态提取基因组路径信息
    genome_paths = extract_genome_paths(state)
    genome_version = state.nextflow_config.get("genome_version")
    star_resource_config = state.resource_config.get("star") if state.resource_config else {}
    
    # 获取来自detect node的基因组信息，为LLM提供判断基础
    genome_info = {}
    if hasattr(state, 'detect_results') and state.detect_results:
        # 从detect结果中获取基因组信息
        detect_results = state.detect_results
        if isinstance(detect_results, dict):
            query_results = detect_results.get('query_results', {})
            if isinstance(query_results, dict) and genome_version:
                # 提取指定基因组的详细信息
                genome_info = query_results.get(genome_version, {})
    
    user_context = {
        "execution_mode": state.execution_mode,
        "genome_config": {
            "genome_version": genome_version,
            "paired_end": state.nextflow_config.get("paired_end")
        },
        "genome_paths": genome_paths,
        "genome_info": genome_info,  # 来自detect node的基因组详细信息
        "star_resource_config": star_resource_config,
        "current_star_params": state.star_params,
        "current_star_index_params": state.star_index_params,
        "fastp_results": state.fastp_results,  # 完整传递FastP结果
        "star_results": state.star_results,  # 历史执行结果
        "optimization_history": state.star_optimization_history
    }

    user_prompt = json.dumps(user_context, ensure_ascii=False, indent=2)
    
    # 创建并调用Agent
    agent_executor = create_star_agent()
    
    # 构建消息格式
    messages = [
        {"role": "user", "content": user_prompt}
    ]
    
    result = await agent_executor.ainvoke({"messages": messages})
    
    # 提取结构化响应
    structured_response = result.get("structured_response")
    try:
        if structured_response:
            log_llm_preview(logger, "star", structured_response)
        else:
            log_llm_preview(logger, "star.raw", {"keys": list(result.keys())[:10]})
    except Exception:
        pass
    
    # 定义最小校验：必须含有 results_dir 和每样本产物文件存在
    def _is_valid_star_results(res: StarResponse) -> bool:
        try:
            sr = getattr(res, 'star_results', {}) or {}
            if not sr.get('success'):
                return False
            results_dir = sr.get('results_dir') or sr.get('results_directory')
            per_outputs = sr.get('per_sample_outputs') or []
            if not results_dir or not per_outputs:
                return False
            missing_paths = []
            for item in per_outputs:
                required = [
                    item.get('aligned_bam'),
                    item.get('log_final'),
                    item.get('log_out'),
                    item.get('log_progress'),
                    item.get('splice_junctions'),
                ]
                if item.get('transcriptome_bam'):
                    required.append(item.get('transcriptome_bam'))
                if item.get('gene_counts'):
                    required.append(item.get('gene_counts'))
                for file_path in required:
                    if file_path and not Path(file_path).exists():
                        missing_paths.append(file_path)
            return len(missing_paths) == 0
        except Exception:
            return False

    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")

    if not _is_valid_star_results(structured_response):
        raise ValueError("Agent返回的结果无效或缺少必要产物")

    return structured_response
