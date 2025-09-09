"""
STAR节点 - 用于执行STAR比对分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, StarResponse
from ..core import get_shared_llm
from ..prompts import STAR_OPTIMIZATION_PROMPT
from ..tools import download_genome_assets, build_star_index, run_nextflow_star, parse_star_metrics, scan_genome_files
import json


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


def star_node(state: AgentState) -> Dict[str, Any]:
    """
    STAR节点实现

    功能：
    - 执行STAR比对
    - 生成比对统计
    - 更新状态信息
    """
    print("\n🎯 STAR比对节点开始执行...")

    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "star" not in completed_steps:
        completed_steps.append("star")

    # 获取执行模式
    execution_mode = state.execution_mode

    # 检查FastP结果依赖（强制）
    if not state.fastp_results or not state.fastp_results.get("success"):
        return {
            "status": "star_failed",
            "response": "❌ STAR比对失败：缺少有效的FastP质控结果，请先完成FastP质量控制",
            "current_step": "star",
            "completed_steps": completed_steps,
        }

    try:
        print(f"⚡ [AGENT] 使用STAR Agent进行比对与资源管理 (模式: {execution_mode})")

        if execution_mode == "single":
            # 单次执行：仅执行比对，不做参数优化
            star_response = _call_star_optimization_agent(state)

            result = {
                "status": "star_completed",
                "response": "✅ STAR比对完成（单次执行模式）\n\n🚀 **执行详情**: 已完成比对，保持原有参数配置",
                "current_step": "star",
                "completed_steps": completed_steps,
                "star_results": {
                    "status": "success",
                    "summary": "STAR比对完成，单次执行模式",
                },
            }
            return result

        elif execution_mode == "optimized":
            # 精细优化：执行+解析+应用优化
            star_response = _call_star_optimization_agent(state)

            result = {
                "status": "star_completed",
                "current_step": "star",
                "completed_steps": completed_steps,
                "star_params": star_response.star_params,
                "star_optimization_suggestions": star_response.star_optimization_suggestions,
                "star_optimization_params": star_response.star_optimization_params,
                "star_results": {
                    "status": "success",
                    "summary": "STAR比对完成，已应用智能优化参数",
                },
            }

            optimization_count = len(star_response.star_optimization_params or {})
            result["response"] = (
                f"✅ STAR比对完成并已优化\n- 比对状态: 成功完成\n- 参数优化: 应用了{optimization_count}个优化参数\n\n"
                f"⚡ **优化详情**: {star_response.star_optimization_suggestions}"
            )
            return result

        elif execution_mode == "batch_optimize":
            # 批次优化：执行+解析+收集优化，不应用
            star_response = _call_star_optimization_agent(state)

            star_optimization = {
                "optimization_reasoning": star_response.star_optimization_suggestions,
                "suggested_params": star_response.star_optimization_params,
                "current_params": state.star_params.copy(),
                "tool_name": "star",
            }

            batch_optimizations = state.batch_optimizations.copy()
            batch_optimizations["star"] = star_optimization

            result = {
                "status": "star_completed",
                "current_step": "star",
                "completed_steps": completed_steps,
                "batch_optimizations": batch_optimizations,
                "star_optimization_suggestions": star_response.star_optimization_suggestions,
                "star_results": {
                    "status": "success",
                    "summary": "STAR比对完成，优化建议已收集",
                },
            }

            optimization_count = len(star_response.star_optimization_params or {})
            result["response"] = (
                f"✅ STAR比对完成\n- 比对状态: 成功完成\n- 优化收集: {optimization_count}个参数优化建议已收集\n\n"
                f"📦 **收集的优化建议**: {star_response.star_optimization_suggestions}"
            )
            return result

        else:
            # 未知模式：按 single 处理
            print(f"ℹ️ 未知执行模式 '{execution_mode}'，按 single 处理")
            star_response = _call_star_optimization_agent(state)
            return {
                "status": "star_completed",
                "response": "✅ STAR比对完成（按single处理）\n\n🚀 **执行详情**: 已完成比对，保持原有参数配置",
                "current_step": "star",
                "completed_steps": completed_steps,
                "star_results": {"status": "success", "summary": "STAR比对完成（single模式）"},
            }

    except Exception as e:
        print(f"❌ STAR节点执行失败: {str(e)}")
        return {
            "status": "star_failed",
            "response": f"❌ STAR比对执行失败: {str(e)}",
            "current_step": "star",
            "completed_steps": completed_steps,
            "star_results": {"status": "failed", "error": str(e)},
        }


def _call_star_optimization_agent(state: AgentState) -> StarResponse:
    """调用STAR Agent，根据执行模式区分行为，统一返回结构化响应"""
    # 组织上下文（仅数据）
    user_context = {
        "execution_mode": state.execution_mode,
        "fastp_results": state.fastp_results,
        "nextflow_config": state.nextflow_config,
        "current_star_params": state.star_params,
        "optimization_history": {
            "star": state.star_optimization_params,
            "fastp": state.fastp_optimization_params,
        },
    }

    # 模式指令
    mode = (state.execution_mode or "single").lower()
    if mode == "single":
        mode_instructions = (
            "本次执行模式为 single（单次执行）。\n"
            "- 仅执行 STAR 比对与必要的资源准备（下载/索引），不进行任何参数优化。\n"
            "- 必须基于 FastP 修剪后的 FASTQ 进行比对。\n"
            "- 保持 current_star_params 原样返回，star_optimization_params 必须为空对象。\n"
            "- 可调用 parse_star_metrics 提取关键指标以形成摘要。\n"
            "- 仍需返回 StarResponse 结构化结果。\n"
        )
    elif mode == "batch_optimize":
        mode_instructions = (
            "本次执行模式为 batch_optimize（批次优化）。\n"
            "- 执行 STAR 并解析结果，生成优化建议，但不要在当前节点应用这些参数。\n"
            "- star_params 请给出“建议后的完整参数字典”，star_optimization_params 仅包含改动的键值对。\n"
        )
    else:  # optimized
        mode_instructions = (
            "本次执行模式为 optimized（精细优化）。\n"
            "- 执行 STAR、解析结果并生成优化建议。\n"
            "- star_params 请返回“应用优化后的完整参数字典”，star_optimization_params 仅包含改动项。\n"
        )

    # 组装用户消息
    user_prompt = (
        "请依据系统提示中的标准流程与指导原则执行本次任务。\n\n"
        + mode_instructions
        + "以下为本次任务的上下文数据（JSON）：\n\n"
        + json.dumps(user_context, ensure_ascii=False, indent=2)
        + "\n\n请基于上述数据完成必要的工具调用，并按系统提示要求返回结构化结果（StarResponse）。"
    )

    # 调用Agent
    agent = create_star_agent()
    messages = [{"role": "user", "content": user_prompt}]
    result = agent.invoke({"messages": messages})

    # 提取结构化响应
    structured = result.get("structured_response") if isinstance(result, dict) else None
    if structured and isinstance(structured, StarResponse):
        return structured

    # 兼容不同返回形态
    if hasattr(result, "content") and isinstance(result.content, StarResponse):
        return result.content
    if hasattr(result, "content") and getattr(result.content, "star_params", None) is not None:
        return result.content

    # 兜底：保持当前参数
    return StarResponse(
        star_params=state.star_params,
        star_optimization_suggestions=(
            "STAR执行完成，单次执行不做优化" if mode == "single" else "STAR执行完成，但未返回有效优化建议"
        ),
        star_optimization_params={},
    )
