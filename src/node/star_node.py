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
    scan_genome_files,
    extract_genome_paths
)
from ..route_decider import decide_next_action_star
from ..logging_bootstrap import get_logger, log_llm_preview
import json
from datetime import datetime

logger = get_logger("rna.nodes.star")


def create_star_agent():
    """创建STAR节点的React Agent"""
    llm = get_shared_llm()

    # 绑定STAR相关工具
    tools = [
        scan_genome_files,
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
        # 统一通过Agent执行STAR，生成优化建议
        logger.info("[STAR] 调用Agent执行STAR比对和优化分析...")
        agent_response = await _call_star_optimization_agent(state)

        # 更新执行参数和优化建议
        optimized_params = agent_response.star_params
        optimization_reasoning = agent_response.star_optimization_suggestions
        optimization_params_changes = agent_response.star_optimization_params

        # 处理执行结果
        star_results = {
            "success": True,
            "status": "success"
        }
        try:
            if getattr(agent_response, 'star_results', None):
                agent_results = agent_response.star_results or {}
                star_results.update(agent_results)
        except Exception:
            star_results["success"] = False
            star_results["status"] = "failed"

        # 生成响应信息
        optimization_count = len(optimization_params_changes or {})
        if optimization_count > 0:
            response = (
                f"✅ STAR比对完成\n- 比对状态: 成功完成\n- 优化分析: 生成了{optimization_count}个优化建议\n\n"
                f"⚡ 优化详情: {optimization_reasoning}"
            )
        else:
            response = (
                "✅ STAR比对完成\n\n"
                "🚀 执行详情: 已完成序列比对，当前参数配置已是最优"
            )
            
        logger.info(f"[STAR] STAR执行完成，生成{optimization_count}个优化参数")

        # 追加优化历史记录
        append_star_optimization_history(
            state=state,
            optimization_params=optimization_params_changes,
            suggestions=optimization_reasoning,
            results=star_results
        )

        # 根据路由决策器结果设置返回上下文
        next_action = decide_next_action_star(state)
        if next_action == "return_confirm":
            state.return_source = "star"
            if not star_results.get("success", True):
                state.return_reason = "failed"
            elif state.execution_mode == 'batch_optimize' and optimization_count > 0:
                state.return_reason = "batch_collect"
            else:
                state.return_reason = "step_confirm"

        # 构建成功结果
        result = {
            "success": True,
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
    
    user_context = {
        "execution_mode": state.execution_mode,
        "genome_config": {
            "genome_version": genome_version,
            "paired_end": state.nextflow_config.get("paired_end")
        },
        "genome_paths": genome_paths,
        "star_resource_config": star_resource_config,
        "current_star_params": state.star_params,
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
    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")
    
    return structured_response
