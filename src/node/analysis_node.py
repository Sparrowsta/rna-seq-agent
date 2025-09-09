"""
Analysis节点 - 用于执行综合分析
"""

from typing import Any, Dict
from langchain_core.messages import HumanMessage
from ..state import AgentState


def analysis_node(state: AgentState) -> Dict[str, Any]:
    """
    Analysis节点实现 - 综合分析节点
    
    功能：
    - 汇总所有分析结果
    - 生成综合报告
    - 提供后续建议
    """
    print("\n📈 综合分析节点开始执行...")
    
    # 获取所有前面步骤的信息
    fastp_results = state.fastp_results or {}
    star_results = state.star_results or {}
    featurecounts_results = state.featurecounts_results or {}
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    species = state.nextflow_config.get('species', 'human')
    
    # 检查是否有足够的数据进行分析
    if not (fastp_results or star_results or featurecounts_results):
        return {
            "status": "error",
            "response": "❌ 缺少分析数据，无法生成综合报告",
        }
    
    # 根据实际执行结果生成综合报告
    analysis_report = f"""
🎉 RNA-seq分析流水线执行完成

📊 **分析概览**:
- 样本数量: {sample_count}个
- 目标物种: {species}
- 流水线: FastP → STAR → FeatureCounts → 分析

✅ **各步骤执行状态**:"""
    
    # 根据实际结果添加各步骤状态
    if fastp_results:
        status_fastp = "✅ 完成" if fastp_results.get("status") == "success" else "❌ 失败"
        analysis_report += f"\n- FastP质控: {status_fastp}"
    
    if star_results:
        status_star = "✅ 完成" if star_results.get("status") == "success" else "❌ 失败"
        analysis_report += f"\n- STAR比对: {status_star}"
        
    if featurecounts_results:
        status_fc = "✅ 完成" if featurecounts_results.get("status") == "success" else "❌ 失败"
        analysis_report += f"\n- 基因定量: {status_fc}"
    
    analysis_report += "\n- 综合分析: ✅ 完成\n"
    
    # 添加分析建议
    analysis_report += """
💡 **后续建议**:
- 可进行差异表达分析
- 建议进行功能富集分析
- 检查结果文件进行进一步分析
    """
    
    print("🧹 [CLEANUP] 清理状态信息，准备下次执行...")
    
    # 返回成功结果并清空影响路由的状态
    return {
        "status": "success",
        "response": analysis_report,
        
        # 清空执行进度状态
        "current_step": "",
        "completed_steps": [],
        "execution_mode": "single",  # 重置为默认模式
        
        # 清空各节点的结果状态
        "fastp_results": {},
        "star_results": {},
        "featurecounts_results": {},
        
        # 清空优化相关状态
        "fastp_optimization_suggestions": "",
        "star_optimization_suggestions": "",
        "featurecounts_optimization_suggestions": "",
    }
