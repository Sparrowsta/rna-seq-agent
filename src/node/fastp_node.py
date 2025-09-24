"""
FastP节点 - 用于执行FastP质量控制
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, FastpResponse
from ..core import get_shared_llm
from ..prompts import FASTP_OPTIMIZATION_PROMPT
from ..tools import (
    run_nextflow_fastp,
    parse_fastp_results
)
from ..route_decider import decide_next_action_fastp
from ..logging_bootstrap import get_logger, log_llm_preview
import json
from datetime import datetime
from pathlib import Path

logger = get_logger("rna.nodes.fastp")


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

def append_fastp_optimization_history(state: AgentState, optimization_params: Dict[str, Any], 
                                    suggestions: str, results: Dict[str, Any]) -> None:
    """追加FastP优化历史记录，保持最近5次记录"""
    history_entry = {
        "timestamp": datetime.now().isoformat(),
        "execution_id": f"fastp_run_{len(state.fastp_optimization_history) + 1}",
        "optimization_params": optimization_params or {},
        "optimization_suggestions": suggestions or "",
        "execution_results": results or {}
    }
    
    # 追加新记录
    state.fastp_optimization_history.append(history_entry)
    
    # 保持最近5次记录
    if len(state.fastp_optimization_history) > 5:
        state.fastp_optimization_history = state.fastp_optimization_history[-5:]
    
    logger.info(f"[FASTP] 已追加优化历史记录，当前保存{len(state.fastp_optimization_history)}次历史")


async def fastp_node(state: AgentState) -> Dict[str, Any]:
    """
    FastP节点实现 - 执行质量控制并生成优化建议
    
    功能：
    - 执行FastP质量控制
    - 基于Agent智能优化参数
    - 更新状态信息
    - 生成optimization_params供路由决策器使用
    """
    logger.info("FastP质控节点开始执行...")
    
    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "fastp" not in completed_steps:
        completed_steps.append("fastp")
    
    try:
        # 统一通过Agent执行FastP质控和优化分析
        logger.info("[FASTP] 调用Agent执行FastP质控和优化分析...")
        agent_response = await _call_fastp_optimization_agent(state)

        # 更新执行参数和优化建议
        optimized_params = agent_response.fastp_params
        optimization_reasoning = agent_response.fastp_optimization_suggestions
        optimization_params_changes = agent_response.fastp_optimization_params

        # 处理执行结果（以工具真实返回为准，默认失败避免"空成功"）
        try:
            agent_results = dict(getattr(agent_response, 'fastp_results', {}) or {})
        except Exception:
            agent_results = {}
        success_flag = bool(agent_results.get("success", False))
        status_text = agent_results.get("status", "success" if success_flag else "failed")
        fastp_results = {**agent_results, "success": success_flag, "status": status_text}

        # 生成响应信息
        optimization_count = len(optimization_params_changes or {})
        if success_flag:
            if optimization_count > 0:
                response = (
                    f"✅ FastP质控完成\n- 质控状态: 成功完成\n- 优化分析: 生成了{optimization_count}个优化建议\n\n"
                    f"⚡ 优化详情: {optimization_reasoning}"
                )
            else:
                response = (
                    "✅ FastP质控完成\n\n"
                    "🚀 执行详情: 已完成质量控制"
                )
        else:
            error_msg = fastp_results.get("error") or fastp_results.get("message") or "FastP执行未产生有效输出"
            response = f"❌ FastP执行失败：{error_msg}"

        logger.info(f"[FASTP] FastP执行完成，生成{optimization_count}个优化参数")

        # 追加优化历史记录
        append_fastp_optimization_history(
            state=state,
            optimization_params=optimization_params_changes,
            suggestions=optimization_reasoning,
            results=fastp_results
        )

        # 更新状态以便路由决策读取最新结果
        state.fastp_results = fastp_results

        # 根据路由决策器结果设置返回上下文
        next_action = decide_next_action_fastp(state)
        if next_action == "return_confirm":
            state.return_source = "fastp"
            if not success_flag:
                state.return_reason = "failed"
            elif state.execution_mode == 'batch_optimize' and optimization_count > 0:
                state.return_reason = "batch_collect"
            else:
                state.return_reason = "step_confirm"

        # 构建成功结果
        result = {
            "success": success_flag,
            "status": "fastp_completed",
            "current_step": "fastp",
            "completed_steps": completed_steps,
            "response": response,
            "fastp_params": optimized_params,
            "fastp_optimization_suggestions": optimization_reasoning,
            "fastp_optimization_params": optimization_params_changes,
            "fastp_results": fastp_results,
        }

        return result

    except Exception as e:
        logger.error(f"[FASTP] FastP执行失败: {str(e)}")

        # 失败时设置返回上下文
        state.return_source = "fastp"
        state.return_reason = "failed"

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


async def _call_fastp_optimization_agent(state: AgentState) -> FastpResponse:
    """调用FastP优化Agent进行智能参数优化"""
    
    # 组织数据上下文（仅数据，不重复流程与指南，遵循系统提示）
    sample_info = {
        "sample_groups": state.nextflow_config.get("sample_groups", []),
        # 补充 paired_end，确保工具拿到与 Prepare 一致的测序类型
        "paired_end": state.nextflow_config.get("paired_end"),
        # 结果目录可选提供，工具内部会自动兜底
        **({"results_dir": state.results_dir} if state.results_dir else {}),
        # M2: 添加 state 信息用于参数版本化
        "state_info": {
            "results_dir": state.results_dir,
            "results_timestamp": state.results_timestamp
        }
    }

    # 提取资源配置片段（仅 FastP）
    fastp_resource_config = state.resource_config.get("fastp") if state.resource_config else {}

    user_context = {
        "execution_mode": state.execution_mode,
        "sample_info": sample_info,
        "fastp_resource_config": fastp_resource_config,
        "current_fastp_params": state.fastp_params,
        "optimization_history": state.fastp_optimization_history
    }

    user_prompt = json.dumps(user_context, ensure_ascii=False, indent=2)
    
    # 创建并调用Agent
    agent_executor = create_fastp_agent()
    
    # 构建消息格式
    messages = [
        {"role": "user", "content": user_prompt}
    ]
    
    result = await agent_executor.ainvoke({"messages": messages})
    
    # 提取结构化响应
    structured_response = result.get("structured_response")
    try:
        if structured_response:
            log_llm_preview(logger, "fastp", structured_response)
        else:
            log_llm_preview(logger, "fastp.raw", {"keys": list(result.keys())[:10]})
    except Exception:
        pass

    # 定义最小校验：必须含有 results_dir 和每样本产物文件存在
    def _is_valid_fastp_results(res: FastpResponse) -> bool:
        try:
            fr = getattr(res, 'fastp_results', {}) or {}
            if not fr.get('success'):
                return False
            results_dir = fr.get('results_dir') or fr.get('results_directory')
            per_outputs = fr.get('per_sample_outputs') or []
            if not results_dir or not per_outputs:
                return False
            missing_paths = []
            for item in per_outputs:
                required = [
                    item.get('html'),
                    item.get('json')
                ]
                # 根据 paired_end 状态添加不同的必需文件
                if item.get('paired_end'):
                    required.extend([
                        item.get('trimmed_r1'),
                        item.get('trimmed_r2')
                    ])
                else:
                    required.append(item.get('trimmed_single'))
                for file_path in required:
                    if file_path and not Path(file_path).exists():
                        missing_paths.append(file_path)
            return len(missing_paths) == 0
        except Exception:
            return False

    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")

    if not _is_valid_fastp_results(structured_response):
        raise ValueError("Agent返回的结果无效或缺少必要产物")

    return structured_response
