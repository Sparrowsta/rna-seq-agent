"""
FeatureCounts节点 - 用于执行基因定量分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FeaturecountsResponse
from ..core import get_shared_llm
from ..prompts import FEATURECOUNTS_OPTIMIZATION_PROMPT
from ..tools import run_nextflow_featurecounts, parse_featurecounts_metrics, scan_genome_files
import json


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


async def featurecounts_node(state: AgentState) -> Dict[str, Any]:
    """
    FeatureCounts节点实现
    
    功能：
    - 执行基因定量
    - 生成表达矩阵  
    - 更新状态信息
    - 根据模式进行参数优化
    """
    print("\n🧬 FeatureCounts定量节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "featurecounts" not in completed_steps:
        completed_steps.append("featurecounts")
    
    # 获取执行模式
    execution_mode = state.execution_mode
    
    if not state.star_results or not state.star_results.get("success"):
        return {
            "success": False,
            "status": "featurecounts_failed",
            "response": "❌ FeatureCounts执行失败：缺少有效的STAR比对结果，请先完成STAR比对",
            "current_step": "featurecounts",
            "completed_steps": completed_steps,
            "featurecounts_results": {
                "success": False,
                "status": "failed",
                "error": "STAR结果不可用或未成功"
            }
        }
    
    try:
        print(f"⚡ [AGENT] 使用FeatureCounts Agent进行定量分析 (模式: {execution_mode})")
        
        if execution_mode == "single":
            # 单次执行：仅执行定量，不做参数优化
            fc_response = await _call_featurecounts_optimization_agent(state)
            
            # 透传Agent返回的results（results_dir, matrix_path等）
            agent_results = getattr(fc_response, 'results', None)
            fc_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                fc_results.update(agent_results)
            
            
            result = {
                "success": True,
                "status": "featurecounts_completed",
                "response": "✅ FeatureCounts定量完成（单次执行模式）\n\n🚀 **执行详情**: 已完成基因定量，保持原有参数配置",
                "current_step": "featurecounts",
                "completed_steps": completed_steps,
                "featurecounts_results": fc_results,
            }
            
            # 同时聚合到跨节点 results 字段，便于统一读取
            try:
                aggregated_results = dict(getattr(state, 'results', {}) or {})
                aggregated_results["featurecounts"] = result.get("featurecounts_results", {})
                result["results"] = aggregated_results
            except Exception:
                pass
                
            return result
            
        elif execution_mode == "optimized":
            # 精细优化：执行+解析+应用优化
            fc_response = await _call_featurecounts_optimization_agent(state)
            
            # 立即更新执行参数
            optimized_params = fc_response.featurecounts_params
            optimization_reasoning = fc_response.featurecounts_optimization_suggestions
            optimization_params_changes = fc_response.featurecounts_optimization_params
            
            # 透传Agent返回的results
            agent_results = getattr(fc_response, 'results', None)
            fc_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                fc_results.update(agent_results)
            
            result = {
                "success": True,
                "status": "featurecounts_completed",
                "current_step": "featurecounts",
                "completed_steps": completed_steps,
                "featurecounts_params": optimized_params,
                "featurecounts_optimization_suggestions": optimization_reasoning,
                "featurecounts_optimization_params": optimization_params_changes,  # 记录变更的参数
                "featurecounts_results": fc_results,
            }
            
            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                f"✅ FeatureCounts定量完成并已优化\n- 定量状态: 成功完成\n- 参数优化: 应用了{optimization_count}个优化参数\n\n"
                f"⚡ **优化详情**: {optimization_reasoning}"
            )
            
            # 同时聚合到跨节点 results 字段，便于统一读取
            try:
                aggregated_results = dict(getattr(state, 'results', {}) or {})
                aggregated_results["featurecounts"] = result.get("featurecounts_results", {})
                result["results"] = aggregated_results
            except Exception:
                pass
                
            return result
            
        elif execution_mode == "yolo":
            # YOLO模式：与optimized相同的执行逻辑，但会自动进入下一步
            fc_response = await _call_featurecounts_optimization_agent(state)
            
            # 立即更新执行参数
            optimized_params = fc_response.featurecounts_params
            optimization_reasoning = fc_response.featurecounts_optimization_suggestions
            optimization_params_changes = fc_response.featurecounts_optimization_params
            
            # 透传Agent返回的results
            agent_results = getattr(fc_response, 'results', None)
            fc_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                fc_results.update(agent_results)
            
            result = {
                "success": True,
                "status": "featurecounts_completed",
                "current_step": "featurecounts",
                "completed_steps": completed_steps,
                "featurecounts_params": optimized_params,
                "featurecounts_optimization_suggestions": optimization_reasoning,
                "featurecounts_optimization_params": optimization_params_changes,
                "featurecounts_results": fc_results,
            }
            
            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                "🎯 FeatureCounts定量完成（YOLO自动模式）\n\n"
                f"⚡ **优化执行**: 已应用{optimization_count}个优化参数，自动进入下一步"
            )
            
            # 同时聚合到跨节点 results 字段，便于统一读取
            try:
                aggregated_results = dict(getattr(state, 'results', {}) or {})
                aggregated_results["featurecounts"] = result.get("featurecounts_results", {})
                result["results"] = aggregated_results
            except Exception:
                pass
                
            return result
            
        elif execution_mode == "batch_optimize":
            # 批次优化：执行+解析+收集优化，不应用
            fc_response = await _call_featurecounts_optimization_agent(state)
            
            # 立即更新参数以供批次收集使用
            optimized_params = fc_response.featurecounts_params
            optimization_reasoning = fc_response.featurecounts_optimization_suggestions
            optimization_params_changes = fc_response.featurecounts_optimization_params
            
            fc_optimization = {
                "optimization_reasoning": optimization_reasoning,
                "suggested_params": optimized_params,
                "optimization_params_changes": optimization_params_changes,
                "current_params": state.featurecounts_params.copy(),
                "tool_name": "featurecounts",
            }
            
            batch_optimizations = state.batch_optimizations.copy()
            batch_optimizations["featurecounts"] = fc_optimization
            
            # 透传Agent返回的results
            agent_results = getattr(fc_response, 'results', None)
            fc_results = {
                "success": True,
                "status": "success",
            }
            if agent_results and isinstance(agent_results, dict):
                fc_results.update(agent_results)
            
            result = {
                "success": True,
                "status": "featurecounts_completed",
                "current_step": "featurecounts",
                "completed_steps": completed_steps,
                "batch_optimizations": batch_optimizations,
                "featurecounts_optimization_suggestions": optimization_reasoning,
                "featurecounts_results": fc_results,
            }
            
            optimization_count = len(optimization_params_changes or {})
            result["response"] = (
                f"✅ FeatureCounts定量完成\n- 定量状态: 成功完成\n- 优化收集: {optimization_count}个参数优化建议已收集\n\n"
                f"📦 **收集的优化建议**: {optimization_reasoning}"
            )
            
            # 同时聚合到跨节点 results 字段，便于统一读取
            try:
                aggregated_results = dict(getattr(state, 'results', {}) or {})
                aggregated_results["featurecounts"] = result.get("featurecounts_results", {})
                result["results"] = aggregated_results
            except Exception:
                pass
                
            return result
            
        else:
            # 未知模式：按 single 处理
            print(f"ℹ️ 未知执行模式 '{execution_mode}'，按 single 处理")
            fc_response = await _call_featurecounts_optimization_agent(state)
            agent_results = getattr(fc_response, 'results', None)
            fc_results = {
                "success": True,
                "status": "success"  # 子结果状态
            }
            if agent_results and isinstance(agent_results, dict):
                fc_results.update(agent_results)
                
            result = {
                "success": True,
                "status": "featurecounts_completed",
                "response": "✅ FeatureCounts定量完成（按single处理）\n\n🚀 **执行详情**: 已完成基因定量，保持原有参数配置",
                "current_step": "featurecounts",
                "completed_steps": completed_steps,
                "featurecounts_results": fc_results,
            }
            
            # 同时聚合到跨节点 results 字段，便于统一读取
            try:
                aggregated_results = dict(getattr(state, 'results', {}) or {})
                aggregated_results["featurecounts"] = result.get("featurecounts_results", {})
                result["results"] = aggregated_results
            except Exception:
                pass
                
            return result
            
    except Exception as e:
        print(f"❌ FeatureCounts节点执行失败: {str(e)}")
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
    """调用FeatureCounts Agent，根据执行模式区分行为，统一返回结构化响应"""
    
    # 组织上下文（仅数据）
    user_context = {
        "execution_mode": state.execution_mode,
        "star_results": state.star_results,
        "nextflow_config": state.nextflow_config,
        "current_featurecounts_params": state.featurecounts_params,
        "genome_version": state.nextflow_config.get("genome_version", ""),
        "optimization_history": {
            "featurecounts": state.featurecounts_optimization_params,
            "star": state.star_optimization_params,
            "fastp": state.fastp_optimization_params,
        },
        # 结果目录相关信息
        **({"results_timestamp": state.results_timestamp} if state.results_timestamp else {}),
        **({"base_results_dir": state.results_dir} if state.results_dir else {}),
    }
    
    # 模式指令
    mode = (state.execution_mode or "single").lower()
    if mode == "single":
        mode_instructions = (
            "本次执行模式为 single（单次执行）。\n"
            "- 仅执行 FeatureCounts 定量，不进行任何参数优化。\n"
            "- 必须基于 STAR 比对后的 BAM 文件进行定量。\n"
            "- 保持 current_featurecounts_params 原样返回（featurecounts_params 可与输入相同），featurecounts_optimization_params 必须为空对象。\n"
            "- 必须调用 run_nextflow_featurecounts 执行，并可调用 parse_featurecounts_metrics 提取关键指标。\n"
            "- 请在结果中返回 results 字段（包含 results_dir、matrix_path 与 per_sample_outputs），便于下游使用。\n"
            "- 仍需返回 FeaturecountsResponse 结构化结果。\n"
        )
    elif mode == "batch_optimize":
        mode_instructions = (
            "本次执行模式为 batch_optimize（批次优化）。\n"
            "- 执行 FeatureCounts 并解析结果，生成优化建议，但不要在当前节点应用这些参数。\n"
            "- featurecounts_params 请给出\"建议后的完整参数字典\"，featurecounts_optimization_params 仅包含改动的键值对。\n"
            "- 返回 results（results_dir, matrix_path, per_sample_outputs）供下游使用。\n"
        )
    else:  # optimized
        mode_instructions = (
            "本次执行模式为 optimized（精细优化）。\n"
            "- 执行 FeatureCounts、解析结果并生成优化建议。\n"
            "- featurecounts_params 请返回\"应用优化后的完整参数字典\"，featurecounts_optimization_params 仅包含改动项。\n"
            "- 返回 results（results_dir, matrix_path, per_sample_outputs）供下游使用。\n"
        )
    
    # 组装用户消息
    user_prompt = (
        "请依据系统提示中的标准流程与指导原则执行本次任务。\n\n"
        + mode_instructions
        + "以下为本次任务的上下文数据（JSON）：\n\n"
        + json.dumps(user_context, ensure_ascii=False, indent=2)
        + "\n\n请基于上述数据完成必要的工具调用，并按系统提示要求返回结构化结果（FeaturecountsResponse）。"
    )
    
    # 创建并调用Agent
    agent_executor = create_featurecounts_agent()
    
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
