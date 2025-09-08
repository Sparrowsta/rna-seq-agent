"""
Analysis节点 - 用于执行综合分析
"""

from typing import Any, Dict
from langchain_core.messages import HumanMessage
from ..state import AgentState


def analysis_node(state: AgentState) -> Dict[str, Any]:
    """
    Analysis节点实现 - 短接版本（测试路由）
    
    功能：
    - 汇总所有分析结果
    - 生成综合报告
    - 提供后续建议
    """
    print("\n📈 综合分析节点开始执行...")
    
    # 短接：直接返回预设信息，不调用LLM
    print("⚡ [SHORT-CIRCUIT] 跳过LLM调用，返回预设结果")
    
    # 获取所有前面步骤的信息
    fastp_results = state.fastp_results or {}
    star_results = state.star_results or {}
    featurecounts_results = state.featurecounts_results or {}
    sample_count = len(state.nextflow_config.get('sample_groups', []))
    species = state.nextflow_config.get('species', 'human')
    
    # 生成综合报告
    analysis_report = f"""
🎉 RNA-seq分析流水线执行完成（短接测试）

📊 **分析概览**:
- 样本数量: {sample_count}个
- 目标物种: {species}
- 流水线: FastP → STAR → FeatureCounts → 分析

✅ **各步骤执行状态**:
- FastP质控: ✅ 通过率95%
- STAR比对: ✅ 比对率88.5%
- 基因定量: ✅ 检测24,587个基因
- 综合分析: ✅ 完成

📋 **主要结果**:
- 高质量读长比例: 92%
- 基因组比对成功率: 88.5%
- 可定量基因数: 24,587个
- 数据质量: 优秀

💡 **后续建议**:
- 可进行差异表达分析
- 建议进行功能富集分析
- 数据质量良好，可用于下游分析
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
        "fastp_optimization_suggestions": "",  # 清空优化建议文字
        "fastp_optimization_params": {},       # 清空优化参数字典
        
        # 可选：保留配置信息供参考，但清空执行状态
        # "nextflow_config": {},  # 保留配置
        # "resource_config": {},   # 保留配置
    }