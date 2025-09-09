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
    
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    execution_mode = state.execution_mode
    
    # 初始化结果字典
    result = {
        "status": "fastp_completed",
        "current_step": "fastp",
        "completed_steps": completed_steps,
        "fastp_results": {
            # 占位，保证后续 update/赋值安全
        }
    }
    
    # 根据执行模式处理优化逻辑
    if execution_mode == "single":
        # 单次执行模式：仅调用 Nextflow FastP，不做结果分析/优化
        print("🚀 [SINGLE] 单次执行模式：调用 Nextflow FastP ...")

        try:
            # 构造样本信息
            sample_groups = state.nextflow_config.get("sample_groups", [])
            if not sample_groups:
                return {
                    "status": "error",
                    "response": "❌ 未找到样本信息(sample_groups)，无法执行FastP",
                    "current_step": "fastp",
                    "completed_steps": completed_steps
                }

            sample_info = {
                "sample_groups": sample_groups,
                **({"results_dir": state.nextflow_config.get("results_dir")} if state.nextflow_config.get("results_dir") else {}),
                **({"results_dir": state.results_dir} if (not state.nextflow_config.get("results_dir") and state.results_dir) else {}),
            }

            fastp_params = state.fastp_params or {}
            fastp_run_result = run_nextflow_fastp.invoke({
                "fastp_params": fastp_params,
                "sample_info": sample_info
            })

            if not (fastp_run_result and fastp_run_result.get("success")):
                error_message = (fastp_run_result or {}).get("error") or "FastP流水线执行失败"
                return {
                    "status": "error",
                    "response": f"❌ {error_message}",
                    "current_step": "fastp",
                    "completed_steps": completed_steps,
                    "fastp_results": {"status": "error", "error": error_message}
                }

            # 成功：填充最小结果集（仅执行信息，不做解析/优化）
            results_dir = fastp_run_result.get("results_dir")
            result["fastp_results"] = {
                "status": "success",
                "results_dir": results_dir,
                "per_sample_outputs": []
            }

            result["response"] = (
                f"✅ FastP质控完成 (单次执行)\n- 样本: {fastp_run_result.get('sample_count')}\n"
                f"- 结果目录: {results_dir}"
            )

        except Exception as e:
            return {
                "status": "error",
                "response": f"❌ FastP单次执行异常: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps
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
            result["fastp_optimization_params"] = optimization_params_changes  # 记录变更的参数

            # 由LLM统一返回的输出路径，直接透传以供下游使用
            try:
                if getattr(agent_response, 'results', None):
                    # 平铺为顶层字段，保持与单次执行一致
                    agent_results = agent_response.results or {}
                    result["fastp_results"] = {
                        "status": "success",
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    }
            except Exception:
                pass

            print(f"✅ [OPTIMIZED] FastP智能优化完成: {len(optimized_params)}个参数")

        except Exception as e:
            print(f"❌ [OPTIMIZED] FastP优化失败: {str(e)}")
            # 直接返回错误，不进行回退处理
            return {
                "status": "error",
                "response": f"❌ FastP智能优化失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps
            }
        
    elif execution_mode == "batch_optimize":
        # 批次优化模式：收集Agent优化参数
        print("📦 [BATCH] FastP批次优化模式，调用Agent收集优化参数...")
        
        try:
            # 调用FastP优化Agent
            agent_response = await _call_fastp_optimization_agent(state)

            # 构建批次优化数据结构
            fastp_optimization = {
                "optimization_reasoning": agent_response.fastp_optimization_suggestions,
                "suggested_params": agent_response.fastp_params,
                "optimization_params_changes": agent_response.fastp_optimization_params,  # 添加变更参数记录
                "current_params": state.fastp_params.copy(),
                "tool_name": "fastp"
            }

            # 将优化参数添加到批次收集器
            batch_optimizations = state.batch_optimizations.copy()
            batch_optimizations["fastp"] = fastp_optimization

            result["batch_optimizations"] = batch_optimizations
            result["response"] += f"\n\n📦 **智能优化参数已收集**: 已收集FastP优化参数"

            # LLM返回的输出路径也一并附带
            try:
                if getattr(agent_response, 'results', None):
                    agent_results = agent_response.results or {}
                    result["fastp_results"] = {
                        "status": "success",
                        "results_dir": agent_results.get("results_dir"),
                        "per_sample_outputs": agent_results.get("per_sample_outputs") or []
                    }
            except Exception:
                pass

            print(f"✅ [BATCH] FastP智能优化参数收集完成: {len(agent_response.fastp_params)}个参数")

        except Exception as e:
            print(f"❌ [BATCH] FastP优化失败: {str(e)}")
            # 直接返回错误，不进行回退处理
            return {
                "status": "error",
                "response": f"❌ FastP批次优化失败: {str(e)}",
                "current_step": "fastp",
                "completed_steps": completed_steps
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

    user_prompt = (
        "请依据系统提示中的标准流程与指导原则进行FastP优化。"
        "以下为本次任务的上下文数据（JSON）：\n\n"
        + json.dumps(user_context, ensure_ascii=False, indent=2)
        + "\n\n请基于上述数据完成工具调用、结果解析与参数优化，并按系统提示要求返回结构化结果。"
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
