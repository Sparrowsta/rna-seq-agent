"""
FastP节点 - 用于执行FastP质量控制
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FastpResponse
from ..core import get_shared_llm
from ..prompts import FASTP_OPTIMIZATION_PROMPT
from ..tools import run_nextflow_fastp, parse_fastp_results
import json


def create_fastp_agent():
    """创建FastP节点的React Agent - 支持真实执行和智能参数优化"""
    llm = get_shared_llm()
    
    # 使用集中管理的系统提示词
    system_prompt = FASTP_OPTIMIZATION_PROMPT
    
    # 集成FastP专用工具
    tools = [
        run_nextflow_fastp,      # 执行FastP质量控制
        parse_fastp_results      # 解析FastP结果文件
    ]
    
    # 创建支持结构化输出的React Agent
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=system_prompt,
        response_format=FastpResponse
    )
    return agent


async def fastp_node(state: AgentState) -> Dict[str, Any]:
    """
    FastP节点实现 - 使用智能Agent进行参数优化
    
    功能：
    - 执行FastP质量控制
    - 基于Agent智能优化参数
    - 更新状态信息
    - 支持批次优化模式
    """
    print("\n🧹 FastP质控节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "fastp" not in completed_steps:
        completed_steps.append("fastp")
    
    execution_mode = state.execution_mode
    
    result = {
        "success": True,
        "status": "fastp_completed",
        "current_step": "fastp",
        "completed_steps": completed_steps,
        "fastp_results": {
            "success": True,
            "status": "success"
        }
    }
    
    if execution_mode == "single":
        # 单次执行：统一通过Agent执行，但不做参数优化
        print("🚀 [SINGLE] 单次执行模式：统一通过Agent执行FastP（不应用优化）")
        try:
            agent_response = await _call_fastp_optimization_agent(state)

            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"].update({
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    })
            except Exception:
                result["success"] = False
                result["status"] = "fastp_failed"  
                result["fastp_results"]["success"] = False
                result["fastp_results"]["status"] = "failed"

            result["response"] = (
                "✅ FastP质控完成（单次执行模式）\n\n"
                "🚀 **执行详情**: 已完成质量控制，保持原有参数配置"
            )
        except Exception as e:
            return {
                "success": False,
                "status": "fastp_failed",
                "response": f"❌ FastP单次执行失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps,
                "fastp_results": {
                    "success": False,
                    "status": "failed", 
                    "error": str(e)
                }
            }
    
    elif execution_mode == "optimized":
        # 精细优化模式：调用Agent进行智能优化
        print("⚡ [OPTIMIZED] 精细优化模式，调用Agent进行智能优化...")
        
        try:
            # 调用FastP优化Agent
            agent_response = await _call_fastp_optimization_agent(state)

            # 立即更新执行参数
            optimized_params = agent_response.fastp_params
            optimization_reasoning = agent_response.fastp_optimization_suggestions
            optimization_params_changes = agent_response.fastp_optimization_params

            result["fastp_params"] = optimized_params
            result["fastp_optimization_suggestions"] = optimization_reasoning
            result["fastp_optimization_params"] = optimization_params_changes

            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"].update({
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    })
            except Exception:
                result["success"] = False
                result["status"] = "fastp_failed"
                result["fastp_results"]["success"] = False
                result["fastp_results"]["status"] = "failed"

            print(f"✅ [OPTIMIZED] FastP智能优化完成: {len(optimized_params)}个参数")

        except Exception as e:
            print(f"❌ [OPTIMIZED] FastP优化失败: {str(e)}")
            return {
                "success": False,
                "status": "fastp_failed",
                "response": f"❌ FastP智能优化失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps,
                "fastp_results": {
                    "success": False,
                    "status": "failed", 
                    "error": str(e)
                }
            }
        
    elif execution_mode == "yolo":
        # YOLO模式：与optimized相同的执行逻辑，但会自动进入下一步
        print("🎯 [YOLO] YOLO模式，自动优化执行...")
        
        try:
            # 调用FastP优化Agent（与optimized相同的逻辑）
            agent_response = await _call_fastp_optimization_agent(state)

            # 立即更新执行参数
            optimized_params = agent_response.fastp_params
            optimization_reasoning = agent_response.fastp_optimization_suggestions
            optimization_params_changes = agent_response.fastp_optimization_params

            result["fastp_params"] = optimized_params
            result["fastp_optimization_suggestions"] = optimization_reasoning
            result["fastp_optimization_params"] = optimization_params_changes

            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"].update({
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    })
            except Exception:
                result["success"] = False
                result["status"] = "fastp_failed"
                result["fastp_results"]["success"] = False
                result["fastp_results"]["status"] = "failed"

            result["response"] = (
                "🎯 FastP质控完成（YOLO自动模式）\n\n"
                "⚡ **优化执行**: 已应用智能参数优化，自动进入下一步"
            )

            print(f"✅ [YOLO] FastP自动优化完成: {len(optimized_params)}个参数")

        except Exception as e:
            print(f"❌ [YOLO] FastP自动优化失败: {str(e)}")
            return {
                "success": False,
                "status": "fastp_failed",
                "response": f"❌ FastP自动优化失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps,
                "fastp_results": {
                    "success": False,
                    "status": "failed", 
                    "error": str(e)
                }
            }
        
    elif execution_mode == "batch_optimize":
        # 批次优化模式：收集Agent优化参数
        print("📦 [BATCH] FastP批次优化模式，调用Agent收集优化参数...")
        
        try:
            # 调用FastP优化Agent
            agent_response = await _call_fastp_optimization_agent(state)

            # 立即更新参数以供批次收集使用
            optimized_params = agent_response.fastp_params
            optimization_reasoning = agent_response.fastp_optimization_suggestions
            optimization_params_changes = agent_response.fastp_optimization_params

            # 构建批次优化数据结构
            fastp_optimization = {
                "optimization_reasoning": optimization_reasoning,
                "suggested_params": optimized_params,
                "optimization_params_changes": optimization_params_changes,  # 添加变更参数记录
                "current_params": state.fastp_params.copy(),
                "tool_name": "fastp"
            }
            # 将优化参数添加到批次收集器
            batch_optimizations = state.batch_optimizations.copy()
            batch_optimizations["fastp"] = fastp_optimization

            result["batch_optimizations"] = batch_optimizations
            result["response"] = (result.get("response", "") + "\n\n📦 **智能优化参数已收集**: 已收集FastP优化参数")

            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"].update({
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    })
            except Exception:
                result["success"] = False
                result["status"] = "fastp_failed"
                result["fastp_results"]["success"] = False
                result["fastp_results"]["status"] = "failed"

            print(f"✅ [BATCH] FastP智能优化参数收集完成: {len(optimized_params)}个参数")

        except Exception as e:
            print(f"❌ [BATCH] FastP优化失败: {str(e)}")
            return {
                "success": False,
                "status": "fastp_failed",
                "response": f"❌ FastP批次优化失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps,
                "fastp_results": {
                    "success": False,
                    "status": "failed",
                    "error": str(e)
                }
            }
    else:
        # 未知模式：按 single 处理
        print(f"ℹ️ 未知执行模式 '{execution_mode}'，按 single 处理")
        try:
            agent_response = await _call_fastp_optimization_agent(state)
            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"].update({
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    })
            except Exception:
                result["success"] = False
                result["status"] = "fastp_failed"
                result["fastp_results"]["success"] = False
                result["fastp_results"]["status"] = "failed"
            result["response"] = (
                "✅ FastP质控完成（按single处理）\n\n"
                "🚀 **执行详情**: 已完成质量控制，保持原有参数配置"
            )
        except Exception as e:
            return {
                "success": False,
                "status": "fastp_failed",
                "response": f"❌ FastP执行失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps,
                "fastp_results": {
                    "success": False,
                    "status": "failed", 
                    "error": str(e)
                }
            }
    
    # 同时聚合到跨节点 results 字段，便于统一读取
    try:
        aggregated_results = dict(getattr(state, 'results', {}) or {})
        aggregated_results["fastp"] = result.get("fastp_results", {})
        result["results"] = aggregated_results
    except Exception:
        pass

    return result


async def _call_fastp_optimization_agent(state: AgentState) -> FastpResponse:
    """调用FastP优化Agent进行智能参数优化"""
    
    # 组织数据上下文（仅数据，不重复流程与指南，遵循系统提示）
    sample_info = {
        "sample_groups": state.nextflow_config.get("sample_groups", []),
        # 结果目录可选提供，工具内部会自动兜底
        **({"results_dir": state.results_dir} if state.results_dir else {})
    }

    user_context = {
        "execution_mode": state.execution_mode,
        "sample_info": sample_info,
        "nextflow_config": state.nextflow_config,
        "current_fastp_params": state.fastp_params,
        "optimization_history": {
            "fastp": state.fastp_optimization_params,
            "star": state.star_optimization_params,
            "featurecounts": state.featurecounts_optimization_params,
        },
    }

    # 根据执行模式，给Agent明确的运行指令，统一走同一个Agent
    mode = (state.execution_mode or "single").lower()
    if mode == "single":
        mode_instructions = (
            "本次执行模式为 single（单次执行）。\n"
            "- 仅执行 FastP 质量控制，不进行任何参数优化。\n"
            "- 保持 current_fastp_params 原样返回（fastp_params 可与输入相同），fastp_optimization_params 必须为空对象。\n"
            "- 必须调用 run_nextflow_fastp 执行，并可调用 parse_fastp_results 解析关键质量指标。\n"
            "- 请在结果中返回 results 字段（包含 results_dir 与 per_sample_outputs），便于下游 STAR 使用。\n"
            "- 仍需返回 FastpResponse 结构化结果。\n"
        )
    elif mode == "batch_optimize":
        mode_instructions = (
            "本次执行模式为 batch_optimize（批次优化）。\n"
            "- 执行 FastP 并解析结果，生成优化建议，但不要在当前节点应用这些参数。\n"
            "- fastp_params 请给出“建议后的完整参数字典”，fastp_optimization_params 仅包含改动的键值对。\n"
            "- 返回 results（results_dir, per_sample_outputs）供下游使用。\n"
        )
    else:  # optimized
        mode_instructions = (
            "本次执行模式为 optimized（精细优化）。\n"
            "- 执行 FastP、解析结果并生成优化建议。\n"
            "- fastp_params 请返回“应用优化后的完整参数字典”，fastp_optimization_params 仅包含改动项。\n"
            "- 返回 results（results_dir, per_sample_outputs）供下游使用。\n"
        )

    user_prompt = (
        "请依据系统提示中的标准流程与指导原则执行本次任务。\n\n"
        + mode_instructions
        + "以下为本次任务的上下文数据（JSON）：\n\n"
        + json.dumps(user_context, ensure_ascii=False, indent=2)
        + "\n\n请基于上述数据完成必要的工具调用，并按系统提示要求返回结构化结果（FastpResponse）。"
    )
    
    # 创建并调用Agent
    agent_executor = create_fastp_agent()
    
    # 构建消息格式
    messages = [
        {"role": "user", "content": user_prompt}
    ]
    
    result = await agent_executor.ainvoke({"messages": messages})
    
    # 提取结构化响应
    structured_response = result.get("structured_response")
    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")
    
    return structured_response
