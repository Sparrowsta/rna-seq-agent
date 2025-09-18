"""
Analysis节点 - LLM驱动的智能分析实现

功能：
- LLM驱动的RNA-seq结果分析
- 智能样本质量评估和健康度判断
- 自动报告生成和文档写入
- 基于create_react_agent的工具调用
"""

from typing import Any, Dict

from langgraph.prebuilt import create_react_agent

from ..state import AgentState
from ..prompts import ANALYSIS_LLM_SYSTEM_PROMPT
from ..tools.analysis_tools import (
    parse_pipeline_results,
    assess_sample_quality,
    write_analysis_report
)
from ..core import get_shared_llm
from ..logging_bootstrap import get_logger

logger = get_logger("rna.analysis")


def analysis_node(state: AgentState) -> Dict[str, Any]:
    """
    Analysis节点 - LLM驱动的智能分析

    使用create_react_agent让LLM自主调用工具：
    1. 解析流水线结果 (parse_pipeline_results)
    2. 评估样本质量 (assess_sample_quality)
    3. 生成分析报告 (write_analysis_report)

    LLM负责核心的分析决策和质量判断
    """
    logger.info("LLM驱动分析节点开始执行")

    try:
        # 构建分析工具列表
        analysis_tools = [
            parse_pipeline_results,
            assess_sample_quality,
            write_analysis_report
        ]

        # 创建LLM Agent
        llm = get_shared_llm()
        agent = create_react_agent(
            llm,
            tools=analysis_tools,
            prompt=ANALYSIS_LLM_SYSTEM_PROMPT
        )

        # 准备分析上下文
        results_dir = _extract_results_directory(state)
        if not results_dir:
            return _create_error_response("无法确定结果目录，分析无法进行")

        sample_groups = state.nextflow_config.get("sample_groups", []) if state.nextflow_config else []
        pipeline_config = state.nextflow_config or {}

        # 构建给LLM的分析请求
        analysis_request = f"""
请对RNA-seq流水线结果进行完整的智能分析：

**任务概述**:
1. 解析结果目录 `{results_dir}` 中的所有流水线步骤结果
2. 评估每个样本的质量和健康度，识别数据模式和问题
3. 生成完整的分析报告并保存到结果目录

**流水线配置**:
- 样本分组: {len(sample_groups)} 个样本配置
- 基因组版本: {pipeline_config.get('genome_version', 'unknown')}
- 物种: {pipeline_config.get('species', 'unknown')}

**你的核心职责**:
1. **数据解析**: 调用 parse_pipeline_results 解析所有步骤结果
2. **质量评估**: 调用 assess_sample_quality 进行智能质量分析和健康度评估
3. **报告生成**: 调用 write_analysis_report 生成最终分析报告

**质量评估标准**:
- 基于RNA-seq专业知识，而非固定阈值
- 考虑数据完整性、一致性和生物学合理性
- 识别异常模式和潜在问题
- 提供可行的优化建议

请开始分析并生成报告。重点关注数据质量洞察和实用建议。
"""

        # 执行LLM Agent分析
        logger.info("启动LLM Agent进行智能分析...")
        messages = [{"role": "user", "content": analysis_request}]

        # Agent执行结果
        agent_result = agent.invoke({"messages": messages})

        # 提取Agent的最终消息作为分析结果
        final_message = agent_result["messages"][-1] if agent_result.get("messages") else ""
        analysis_response = ""
        if hasattr(final_message, 'content'):
            analysis_response = final_message.content
        else:
            analysis_response = str(final_message)

        logger.info("LLM分析完成，生成用户响应")

        # 构建成功响应
        user_response = f"""
🎉 RNA-seq智能分析完成！

🤖 AI分析师工作总结:
{analysis_response}

📁 分析结果已保存到: {results_dir}
💡 请查看生成的分析报告了解详细指标和建议
"""

        return {
            "success": True,
            "status": "analysis_completed",
            "response": user_response.strip(),
            "analysis_agent_result": agent_result,
            "analysis_report_path": results_dir,

            # 清空执行进度状态
            "current_step": "",
            "completed_steps": [],
            "execution_mode": "single",

            # 清空各节点的结果状态
            "fastp_results": {},
            "star_results": {},
            "hisat2_results": {},
            "featurecounts_results": {},

            # 标记RNA-seq流程完成
            "rna_seq_complete": True
        }

    except Exception as e:
        logger.error(f"LLM分析节点执行失败: {str(e)}", exc_info=True)
        return _create_error_response(f"LLM分析节点执行失败: {str(e)}")


def _extract_results_directory(state: AgentState) -> str:
    """从状态中提取结果目录路径"""
    # 优先从query_results中获取
    if state.query_results:
        return state.query_results.get("results_dir", "")

    # 回退到results_dir字段
    if hasattr(state, 'results_dir') and state.results_dir:
        return state.results_dir

    # 最后从各执行结果中推断
    for result_key in ['fastp_results', 'star_results', 'hisat2_results', 'featurecounts_results']:
        result = getattr(state, result_key, {})
        if result and result.get("results_directory"):
            return result["results_directory"]

    return ""


def _create_error_response(error_message: str) -> Dict[str, Any]:
    """创建错误响应"""
    return {
        "success": False,
        "status": "analysis_error",
        "response": f"❌ 分析失败: {error_message}",
        "analysis_report_path": "",
        "rna_seq_complete": False
    }