"""
HISAT2节点 - 用于执行HISAT2比对分析
"""

from typing import Any, Dict
from langgraph.prebuilt import create_react_agent
from ..state import AgentState, Hisat2Response
from ..core import get_shared_llm
from ..prompts import HISAT2_OPTIMIZATION_PROMPT
from ..tools import (
    download_genome_assets,
    scan_genome_files,
    build_hisat2_index,
    run_nextflow_hisat2,
    parse_hisat2_metrics
)
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


def append_hisat2_optimization_history(state: AgentState, optimization_params: Dict[str, Any]) -> None:
    """追加HISAT2优化历史记录"""
    if not hasattr(state, 'hisat2_optimization_history') or state.hisat2_optimization_history is None:
        state.hisat2_optimization_history = []

    # 追加当前优化记录
    state.hisat2_optimization_history.append(optimization_params)

    # 限制历史记录数量，只保留最近5次
    if len(state.hisat2_optimization_history) > 5:
        state.hisat2_optimization_history = state.hisat2_optimization_history[-5:]

    logger.info(f"[HISAT2] 已追加优化历史记录，当前保存{len(state.hisat2_optimization_history)}次历史")
    logger.info(f"[HISAT2] 最新优化历史记录：{state.hisat2_optimization_history[-1]}")


async def hisat2_node(state: AgentState) -> Dict[str, Any]:
    """
    HISAT2节点实现 - 执行序列比对并生成优化建议
    
    功能：
    - 执行HISAT2比对
    - 基于Agent智能优化参数
    - 更新状态信息
    - 生成optimization_params供路由决策器使用
    """
    logger.info("HISAT2比对节点开始执行...")

    # 更新执行进度
    completed_steps = state.completed_steps.copy() if state.completed_steps else []
    if "hisat2" not in completed_steps:
        completed_steps.append("hisat2")

    # 依赖检查：需要FastP成功结果
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
        # 统一通过Agent执行HISAT2，生成优化建议
        logger.info("[HISAT2] 调用Agent执行HISAT2比对和优化分析...")
        agent_response = await _call_hisat2_optimization_agent(state)

        # 更新执行参数和优化建议
        optimized_params = agent_response.hisat2_params
        optimization_reasoning = agent_response.hisat2_optimization_suggestions
        optimization_params_changes = agent_response.hisat2_optimization_params

        # 处理执行结果
        hisat2_results = {
            "success": True,
            "status": "success"
        }
        try:
            if getattr(agent_response, 'results', None):
                agent_results = agent_response.results or {}
                hisat2_results.update(agent_results)
        except Exception:
            hisat2_results["success"] = False
            hisat2_results["status"] = "failed"

        # 生成响应信息
        optimization_count = len(optimization_params_changes or {})
        if optimization_count > 0:
            response = (
                f"✅ HISAT2比对完成\n- 比对状态: 成功完成\n- 优化分析: 生成了{optimization_count}个优化建议\n\n"
                f"⚡ 优化详情: {optimization_reasoning}"
            )
        else:
            response = (
                "✅ HISAT2比对完成\n\n"
                "🚀 执行详情: 已完成序列比对，当前参数配置已是最优"
            )

        logger.info(f"[HISAT2] HISAT2执行完成，生成{optimization_count}个优化参数")

        # 构建成功结果
        result = {
            "success": True,
            "status": "hisat2_completed",
            "current_step": "hisat2",
            "completed_steps": completed_steps,
            "response": response,
            "hisat2_params": optimized_params,
            "hisat2_optimization_suggestions": optimization_reasoning,
            "hisat2_optimization_params": optimization_params_changes,
            "hisat2_results": hisat2_results,
        }



        return result

    except Exception as e:
        logger.error(f"[HISAT2] HISAT2执行失败: {str(e)}")
        return {
            "success": False,
            "status": "hisat2_failed",
            "response": f"❌ HISAT2执行失败: {str(e)}",
            "current_step": "hisat2",
            "completed_steps": completed_steps,
            "hisat2_results": {
                "success": False,
                "status": "failed", 
                "error": str(e)
            }
        }


async def _call_hisat2_optimization_agent(state: AgentState) -> Hisat2Response:
    """调用HISAT2优化Agent进行智能参数优化"""

    # 获取基因组配置信息（从detect节点的query_results中获取）
    genome_version = state.nextflow_config.get("genome_version")
    genome_info = {}
    if genome_version and state.query_results:
        try:
            # 从detect节点已经收集的结果中获取基因组配置
            genome_scan_result = state.query_results.get("verify_genome_setup", {})
            if isinstance(genome_scan_result, dict):
                genomes_dict = genome_scan_result.get("genomes", {})
                if isinstance(genomes_dict, dict) and genome_version in genomes_dict:
                    genome_info = genomes_dict[genome_version]
        except Exception as e:
            logger.warning(f"从query_results获取基因组配置失败: {e}")

    # 组织数据上下文（仅数据，不重复流程与指南，遵循系统提示）
    sample_info = {
        "sample_groups": state.nextflow_config.get("sample_groups", []),
        # 结果目录可选提供，工具内部会自动兜底
        **({"results_dir": state.results_dir} if state.results_dir else {}),
        # 添加state信息用于参数版本化
        "state_info": {
            "results_dir": state.results_dir,
            "results_timestamp": state.results_timestamp
        }
    }

    user_context = {
        "execution_mode": state.execution_mode,
        "genome_config": {
            "genome_version": genome_version,
            "paired_end": state.nextflow_config.get("paired_end")
        },
        "genome_info": genome_info,
        "sample_info": sample_info,
        "current_hisat2_params": state.hisat2_params,
        "fastp_results": state.fastp_results,
        "hisat2_results": state.hisat2_results,  # 新增：历史执行结果
        "optimization_history": {
            "hisat2": state.hisat2_optimization_history,  # 完整历史列表
            "fastp": state.fastp_optimization_params,
        },
    }

    user_prompt = json.dumps(user_context, ensure_ascii=False, indent=2)
    
    # 创建并调用Agent
    agent_executor = create_hisat2_agent()
    
    # 构建消息格式
    messages = [
        {"role": "user", "content": user_prompt}
    ]
    
    result = await agent_executor.ainvoke({"messages": messages})
    
    # 提取结构化响应
    structured_response = result.get("structured_response")
    try:
        if structured_response:
            log_llm_preview(logger, "hisat2", structured_response)
        else:
            log_llm_preview(logger, "hisat2.raw", {"keys": list(result.keys())[:10]})
    except Exception:
        pass
    if not structured_response:
        raise ValueError("Agent返回的结构化响应为空")
    
    return structured_response
