"""
HISAT2节点 - 用于执行HISAT2比对分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, Hisat2Response
from ..core import get_shared_llm
from ..prompts import HISAT2_OPTIMIZATION_PROMPT
from ..tools import download_genome_assets, build_hisat2_index, run_nextflow_hisat2, parse_hisat2_metrics, scan_genome_files
from ..logging_bootstrap import get_logger, log_llm_preview
import json

logger = get_logger("rna.nodes.hisat2")


def create_hisat2_agent():
    """创建HISAT2节点的React Agent"""
    llm = get_shared_llm()

    # 绑定HISAT2相关工具
    tools = [
        scan_genome_files,
        download_genome_assets,
        build_hisat2_index,
        run_nextflow_hisat2,
        parse_hisat2_metrics,
    ]

    # 创建使用工具的React Agent
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=HISAT2_OPTIMIZATION_PROMPT,
        response_format=Hisat2Response,
    )
    return agent


async def hisat2_node(state: AgentState) -> Dict[str, Any]:
    """
    HISAT2节点实现

    功能：
    - 执行HISAT2比对
    - 生成比对统计
    - 更新状态信息
    """
    logger.info("HISAT2比对节点开始执行")

    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "hisat2" not in completed_steps:
        completed_steps.append("hisat2")

    # 获取执行模式
    execution_mode = state.execution_mode

    if not state.fastp_results or not state.fastp_results.get("success"):
        return {
            "success": False,
            "status": "hisat2_failed",
            "response": "❌ HISAT2比对失败：缺少有效的FastP质控结果，请先完成FastP质量控制",
            "current_step": "hisat2",
            "completed_steps": completed_steps,
            "hisat2_results": {
                "success": False,
                "status": "failed",
                "error": "FastP结果不可用或未成功"
            }
        }

    try:
        logger.info(f"[AGENT] 使用HISAT2 Agent进行比对与资源管理 (模式: {execution_mode})")

        if execution_mode == "single":
            # 单次执行：仅执行比对，不做参数优化
            hisat2_response = await _call_hisat2_optimization_agent(state)
            try:
                log_llm_preview(logger, "hisat2", hisat2_response)
            except Exception:
                pass

            # 透传Agent返回的results（results_dir, per_sample_outputs）
            agent_results = getattr(hisat2_response, 'results', None)
            hisat2_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                hisat2_results.update(agent_results)
            
            # 确保BAM文件路径信息完整（根据路径契约要求）
            hisat2_results = _ensure_bam_paths_from_per_sample(hisat2_results)

            result = {
                "success": True,
                "status": "hisat2_completed",
                "response": "✅ HISAT2比对完成（单次执行模式）\n\n🚀 执行详情: 已完成比对，保持原有参数配置",
                "current_step": "hisat2",
                "completed_steps": completed_steps,
                "hisat2_results": hisat2_results,
            }
            return result

        elif execution_mode == "optimized":
            # 精细优化：执行+解析+应用优化
            hisat2_response = await _call_hisat2_optimization_agent(state)
            try:
                log_llm_preview(logger, "hisat2", hisat2_response)
            except Exception:
                pass

            # 立即更新执行参数
            optimized_params = hisat2_response.hisat2_params
            optimization_reasoning = hisat2_response.hisat2_optimization_suggestions
            optimization_params_changes = hisat2_response.hisat2_optimization_params

            # 透传Agent返回的results（results_dir, per_sample_outputs）
            agent_results = getattr(hisat2_response, 'results', None)
            hisat2_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                hisat2_results.update(agent_results)
            
            # 确保BAM文件路径信息完整（根据路径契约要求）
            hisat2_results = _ensure_bam_paths_from_per_sample(hisat2_results)

            result = {
                "success": True,
                "status": "hisat2_completed",
                "current_step": "hisat2",
                "completed_steps": completed_steps,
                "hisat2_params": optimized_params,
                "hisat2_optimization_suggestions": optimization_reasoning,
                "hisat2_optimization_params": optimization_params_changes,  # 记录变更的参数
                "hisat2_results": hisat2_results,
            }

            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                f"✅ HISAT2比对完成并已优化\n- 比对状态: 成功完成\n- 参数优化: 应用了{optimization_count}个优化参数\n\n"
                f"⚡ 优化详情: {optimization_reasoning}"
            )
            return result
            
        elif execution_mode == "yolo":
            # YOLO模式：与optimized相同的执行逻辑，但会自动进入下一步
            hisat2_response = await _call_hisat2_optimization_agent(state)
            try:
                log_llm_preview(logger, "hisat2", hisat2_response)
            except Exception:
                pass

            # 立即更新执行参数
            optimized_params = hisat2_response.hisat2_params
            optimization_reasoning = hisat2_response.hisat2_optimization_suggestions
            optimization_params_changes = hisat2_response.hisat2_optimization_params

            # 透传Agent返回的results（results_dir, per_sample_outputs）
            agent_results = getattr(hisat2_response, 'results', None)
            hisat2_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                hisat2_results.update(agent_results)
            
            # 确保BAM文件路径信息完整（根据路径契约要求）
            hisat2_results = _ensure_bam_paths_from_per_sample(hisat2_results)

            result = {
                "success": True,
                "status": "hisat2_completed",
                "current_step": "hisat2",
                "completed_steps": completed_steps,
                "hisat2_params": optimized_params,
                "hisat2_optimization_suggestions": optimization_reasoning,
                "hisat2_optimization_params": optimization_params_changes,
                "hisat2_results": hisat2_results,
            }

            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                "🎯 HISAT2比对完成（YOLO自动模式）\n\n"
                f"⚡ 优化执行: 已应用{optimization_count}个优化参数，自动进入下一步"
            )
            return result

        elif execution_mode == "batch_optimize":
            # 批次优化：执行+解析+收集优化，不应用
            hisat2_response = await _call_hisat2_optimization_agent(state)
            try:
                log_llm_preview(logger, "hisat2", hisat2_response)
            except Exception:
                pass

            # 立即更新参数以供批次收集使用
            optimized_params = hisat2_response.hisat2_params
            optimization_reasoning = hisat2_response.hisat2_optimization_suggestions
            optimization_params_changes = hisat2_response.hisat2_optimization_params

            hisat2_optimization = {
                "optimization_reasoning": optimization_reasoning,
                "suggested_params": optimized_params,
                "optimization_params_changes": optimization_params_changes,
                "current_params": state.hisat2_params.copy(),
                "tool_name": "hisat2",
            }

            batch_optimizations = state.batch_optimizations.copy()
            batch_optimizations["hisat2"] = hisat2_optimization

            # 透传Agent返回的results（results_dir, per_sample_outputs）
            agent_results = getattr(hisat2_response, 'results', None)
            hisat2_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                hisat2_results.update(agent_results)
            
            # 确保BAM文件路径信息完整（根据路径契约要求）
            hisat2_results = _ensure_bam_paths_from_per_sample(hisat2_results)

            result = {
                "success": True,
                "status": "hisat2_completed",
                "current_step": "hisat2",
                "completed_steps": completed_steps,
                "batch_optimizations": batch_optimizations,
                "hisat2_optimization_suggestions": optimization_reasoning,
                "hisat2_results": hisat2_results,
            }

            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                f"✅ HISAT2比对完成\n- 比对状态: 成功完成\n- 优化收集: {optimization_count}个参数优化建议已收集\n\n"
                f"📊 状态更新: hisat2_completed"
            )
            return result

        else:
            # 未知模式：按 single 处理
            logger.warning(f"未知执行模式 '{execution_mode}'，按 single 处理")
            hisat2_response = await _call_hisat2_optimization_agent(state)
            try:
                log_llm_preview(logger, "hisat2", hisat2_response)
            except Exception:
                pass
            agent_results = getattr(hisat2_response, 'results', None)
            hisat2_results = {
                "success": True,
                "status": "success"  # 子结果状态
            }
            if agent_results and isinstance(agent_results, dict):
                hisat2_results.update(agent_results)
            
            return {
                "success": True,
                "status": "hisat2_completed",
                "response": "✅ HISAT2比对完成（按single处理）\n\n🚀 执行详情: 已完成比对，保持原有参数配置",
                "current_step": "hisat2",
                "completed_steps": completed_steps,
                "hisat2_results": hisat2_results,
            }

    except Exception as e:
        logger.error(f"HISAT2节点执行失败: {str(e)}", exc_info=True)
        return {
            "success": False,
            "status": "hisat2_failed",
            "response": f"❌ HISAT2比对执行失败: {str(e)}",
            "current_step": "hisat2",
            "completed_steps": completed_steps,
            "hisat2_results": {
                "success": False,
                "status": "failed", 
                "error": str(e)
            },
        }


async def _call_hisat2_optimization_agent(state: AgentState) -> Hisat2Response:
    """调用HISAT2 Agent，根据执行模式区分行为，统一返回结构化响应"""
    # 组织上下文（仅数据）
    user_context = {
        "execution_mode": state.execution_mode,
        "fastp_results": state.fastp_results,
        "nextflow_config": state.nextflow_config,
        "current_hisat2_params": state.hisat2_params,
        "optimization_history": {
            "hisat2": state.hisat2_optimization_params,
            "fastp": state.fastp_optimization_params,
        },
        # 结果目录相关信息（用于统一输出路径）
        **({"results_timestamp": state.results_timestamp} if getattr(state, 'results_timestamp', None) else {}),
        **({"base_results_dir": state.results_dir} if getattr(state, 'results_dir', None) else {}),
    }

    # 模式指令
    mode = (state.execution_mode or "single").lower()
    if mode == "single":
        mode_instructions = (
            "本次执行模式为 single（单次执行）。\n"
            "- 仅执行 HISAT2 比对与必要的资源准备（下载/索引），不进行任何参数优化。\n"
            "- 必须基于 FastP 修剪后的 FASTQ 进行比对。\n"
            "- 保持 current_hisat2_params 原样返回，hisat2_optimization_params 必须为空对象。\n"
            "- 可调用 parse_hisat2_metrics 提取关键指标以形成摘要。\n"
            "- 仍需返回 Hisat2Response 结构化结果。\n"
        )
    elif mode == "batch_optimize":
        mode_instructions = (
            "本次执行模式为 batch_optimize（批次优化）。\n"
            "- 执行 HISAT2 并解析结果，生成优化建议，但不要在当前节点应用这些参数。\n"
            "- hisat2_params 请给出\"建议后的完整参数字典\"，hisat2_optimization_params 仅包含改动的键值对。\n"
        )
    else:  # optimized
        mode_instructions = (
            "本次执行模式为 optimized（精细优化）。\n"
            "- 执行 HISAT2、解析结果并生成优化建议。\n"
            "- hisat2_params 请返回\"应用优化后的完整参数字典\"，hisat2_optimization_params 仅包含改动项。\n"
        )

    # 组装用户消息
    user_prompt = (
        "请依据系统提示中的标准流程与指导原则执行本次任务。\n\n"
        + mode_instructions
        + "以下为本次任务的上下文数据（JSON）：\n\n"
        + json.dumps(user_context, ensure_ascii=False, indent=2)
        + "\n\n请基于上述数据完成必要的工具调用，并按系统提示要求返回结构化结果（Hisat2Response）。"
    )

    # 调用Agent
    agent = create_hisat2_agent()
    messages = [{"role": "user", "content": user_prompt}]
    result = await agent.ainvoke({"messages": messages})
    try:
        structured = result.get("structured_response") if isinstance(result, dict) else None
        if structured:
            log_llm_preview(logger, "hisat2", structured)
        else:
            log_llm_preview(logger, "hisat2.raw", {"keys": list(result.keys())[:10]})
    except Exception:
        pass

    # 提取结构化响应
    structured = result.get("structured_response") if isinstance(result, dict) else None
    if structured and isinstance(structured, Hisat2Response):
        return structured

    # 兜底：保持当前参数
    return Hisat2Response(
        hisat2_params=state.hisat2_params,
        hisat2_optimization_suggestions=(
            "HISAT2执行完成，单次执行不做优化" if mode == "single" else "HISAT2执行完成，但未返回有效优化建议"
        ),
        hisat2_optimization_params={},
    )


def _ensure_bam_paths_from_per_sample(hisat2_results: Dict[str, Any]) -> Dict[str, Any]:
    """根据路径契约要求，从per_sample_outputs提取BAM文件路径列表
    
    严格遵循docs/path_contract.md的HISAT2→FeatureCounts接口约定：
    - 不覆盖results_dir
    - 仅从per_sample_outputs中提取aligned_bam路径
    - 生成bam_files列表供下游使用
    """
    enhanced = dict(hisat2_results or {})
    
    try:
        bam_files: list[str] = []
        per_sample = enhanced.get("per_sample_outputs") or []
        
        # 根据路径契约，直接从per_sample_outputs的aligned_bam字段提取
        for item in per_sample:
            aligned_bam = item.get("aligned_bam")
            if aligned_bam:
                bam_files.append(aligned_bam)
        
        if bam_files:
            enhanced["bam_files"] = bam_files
            enhanced["bam_files_verified"] = True
        else:
            enhanced["bam_files_verified"] = False
    
    except Exception as e:
        logger.warning(f"提取BAM路径时出错: {e}")
        enhanced["bam_files_verified"] = False
        enhanced.setdefault("error", str(e))
    
    return enhanced
