import json
import re
from typing import Dict, Any
from pathlib import Path

from ..core import get_shared_llm
from ..state import AgentState, AnalysisResponse


def extract_nextflow_metrics(data_dir: str = ".") -> Dict[str, Any]:
    """提取Nextflow执行指标 - 直接读取报告文件"""
    metrics = {}
    
    # 读取执行报告
    report_file = Path(data_dir) / "reports" / "execution_report.txt"
    timeline_file = Path(data_dir) / "reports" / "execution_timeline.txt"
    trace_file = Path(data_dir) / "reports" / "execution_trace.txt"
    
    try:
        # 1. 读取执行报告
        if report_file.exists():
            with open(report_file, 'r', encoding='utf-8') as f:
                report_content = f.read()
                
            # 提取工作流状态信息
            if "Workflow completed successfully" in report_content:
                metrics["workflow_status"] = "success"
            elif "Workflow execution failed" in report_content or "ERROR" in report_content:
                metrics["workflow_status"] = "failed"
            else:
                metrics["workflow_status"] = "unknown"
                
        # 2. 读取执行跟踪
        if trace_file.exists():
            with open(trace_file, 'r', encoding='utf-8') as f:
                trace_content = f.read()
                
            # 统计进程执行情况
            lines = trace_content.strip().split('\n')
            if len(lines) > 1:  # 跳过表头
                process_stats = {"total": 0, "completed": 0, "failed": 0}
                for line in lines[1:]:
                    parts = line.split('\t')
                    if len(parts) > 6:  # 确保有足够的列
                        status = parts[6].strip()  # status列
                        process_stats["total"] += 1
                        if status == "COMPLETED":
                            process_stats["completed"] += 1
                        elif status == "FAILED":
                            process_stats["failed"] += 1
                            
                metrics["process_stats"] = process_stats
                
    except (FileNotFoundError, UnicodeDecodeError, IndexError) as e:
        # 如果报告文件不存在或读取失败，回退到从执行输出推断
        pass
    
    return metrics


def extract_fastp_metrics(data_dir: str = ".") -> Dict[str, Any]:
    """提取fastp质控指标 - 直接读取JSON报告文件"""
    metrics = {}
    
    fastp_dir = Path(data_dir) / "results" / "fastp"
    if not fastp_dir.exists():
        return metrics
    
    # 合并所有样本的fastp统计
    total_before_reads = 0
    total_after_reads = 0
    total_before_bases = 0
    total_after_bases = 0
    sample_count = 0
    
    for sample_dir in fastp_dir.iterdir():
        if not sample_dir.is_dir():
            continue
            
        json_files = list(sample_dir.glob("*.fastp.json"))
        if json_files:
            try:
                with open(json_files[0], 'r') as f:
                    data = json.load(f)
                    summary = data.get("summary", {})
                    
                    before = summary.get("before_filtering", {})
                    after = summary.get("after_filtering", {})
                    
                    total_before_reads += before.get("total_reads", 0)
                    total_after_reads += after.get("total_reads", 0)
                    total_before_bases += before.get("total_bases", 0)
                    total_after_bases += after.get("total_bases", 0)
                    sample_count += 1
                    
            except (json.JSONDecodeError, FileNotFoundError, KeyError):
                continue
    
    if sample_count > 0:
        metrics["total_samples"] = sample_count
        metrics["total_before_reads"] = total_before_reads
        metrics["total_after_reads"] = total_after_reads
        metrics["total_before_bases"] = total_before_bases  
        metrics["total_after_bases"] = total_after_bases
        metrics["filtering_rate"] = f"{((total_before_reads - total_after_reads) / total_before_reads * 100):.1f}%" if total_before_reads > 0 else "0%"
    
    return metrics


def extract_star_metrics(execution_output: str) -> Dict[str, Any]:
    """提取STAR比对指标"""
    metrics = {}
    
    patterns = {
        "uniquely_mapped": r'Uniquely mapped reads %\s+\|\s+(\d+\.?\d*)%',
        "multimapped": r'% of reads mapped to multiple loci\s+\|\s+(\d+\.?\d*)%',
        "unmapped": r'% of reads unmapped.*?\s+\|\s+(\d+\.?\d*)%'
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, execution_output)
        if match:
            metrics[key] = match.group(1) + "%"
    
    return metrics


def extract_featurecounts_metrics(data_dir: str = ".") -> Dict[str, Any]:
    """提取featureCounts定量指标 - 直接读取summary文件"""
    metrics = {}
    
    summary_file = Path(data_dir) / "results" / "featurecounts" / "all_samples.counts.txt.summary"
    if not summary_file.exists():
        return metrics
    
    try:
        with open(summary_file, 'r') as f:
            lines = f.readlines()
            
        if len(lines) < 2:
            return metrics
            
        # 第一行是表头，包含样本名
        header = lines[0].strip().split('\t')
        sample_count = len(header) - 1  # 减去Status列
        
        total_assigned = 0
        total_unassigned = 0
        
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
                
            status = parts[0]
            counts = [int(x) for x in parts[1:] if x.isdigit()]
            
            if status == "Assigned":
                total_assigned = sum(counts)
            elif status.startswith("Unassigned"):
                total_unassigned += sum(counts)
        
        total_reads = total_assigned + total_unassigned
        
        metrics["sample_count"] = sample_count
        metrics["assigned_reads"] = total_assigned
        metrics["unassigned_reads"] = total_unassigned
        metrics["total_processed_reads"] = total_reads
        metrics["assignment_rate"] = f"{(total_assigned / total_reads * 100):.1f}%" if total_reads > 0 else "0%"
        
    except (FileNotFoundError, ValueError, IndexError):
        pass
    
    return metrics


def save_analysis_report(analysis_response: AnalysisResponse, data_dir: str = ".") -> Dict[str, str]:
    """保存分析报告到文件"""
    import datetime
    
    results_path = Path(data_dir) / "results"
    results_path.mkdir(exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    saved_files = {}
    
    try:
        # 1. 保存JSON格式的完整数据
        json_file = results_path / f"analysis_report_{timestamp}.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump({
                "timestamp": timestamp,
                "analysis_summary": analysis_response.analysis_summary,
                "analysis_insights": analysis_response.analysis_insights,
                "result_files": analysis_response.result_files,
                "quality_metrics": analysis_response.quality_metrics,
                "next_steps": analysis_response.next_steps
            }, f, ensure_ascii=False, indent=2)
        saved_files["json"] = str(json_file.relative_to(Path(data_dir)))
        
        # 2. 保存Markdown格式的可读报告
        md_file = results_path / f"analysis_summary_{timestamp}.md"
        with open(md_file, 'w', encoding='utf-8') as f:
            f.write(f"# RNA-seq 分析报告\n\n")
            f.write(f"**生成时间**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write(f"## 分析总结\n\n{analysis_response.analysis_summary}\n\n")
            
            if analysis_response.analysis_insights:
                f.write(f"## 分析洞察\n\n")
                for insight in analysis_response.analysis_insights:
                    f.write(f"- {insight}\n")
                f.write(f"\n")
            
            if analysis_response.result_files:
                f.write(f"## 结果文件\n\n")
                for category, path in analysis_response.result_files.items():
                    f.write(f"- **{category}**: `{path}`\n")
                f.write(f"\n")
            
            if analysis_response.next_steps:
                f.write(f"## 下一步建议\n\n")
                for step in analysis_response.next_steps:
                    f.write(f"- {step}\n")
        saved_files["markdown"] = str(md_file.relative_to(Path(data_dir)))
        
        # 3. 保存最新报告的符号链接
        latest_json = results_path / "analysis_report_latest.json"
        latest_md = results_path / "analysis_summary_latest.md"
        
        # 删除旧的符号链接(如果存在)
        if latest_json.is_symlink():
            latest_json.unlink()
        if latest_md.is_symlink():
            latest_md.unlink()
            
        # 创建新的符号链接
        latest_json.symlink_to(json_file.name)
        latest_md.symlink_to(md_file.name)
        
        saved_files["latest_json"] = str(latest_json.relative_to(Path(data_dir)))
        saved_files["latest_md"] = str(latest_md.relative_to(Path(data_dir)))
        
    except Exception as e:
        print(f"保存分析报告时出错: {e}")
    
    return saved_files


def scan_result_files(data_dir: str = ".") -> Dict[str, Any]:
    """扫描并统计结果文件"""
    file_stats = {
        "total_files": 0,
        "file_categories": {},
        "key_files": {}
    }
    
    results_path = Path(data_dir) / "results"
    if not results_path.exists():
        return file_stats
    
    # 定义文件类别
    categories = {
        "quality_reports": [".html", ".json"],
        "alignment_files": [".bam", ".sam"],
        "quantification": [".txt", ".tsv", ".counts"],
        "logs": [".log", ".out", ".err"]
    }
    
    all_files = list(results_path.rglob("*"))
    file_stats["total_files"] = len([f for f in all_files if f.is_file()])
    
    for category, extensions in categories.items():
        category_files = []
        for ext in extensions:
            category_files.extend([f for f in all_files if f.suffix == ext])
        
        file_stats["file_categories"][category] = len(category_files)
        if category_files:
            # 记录第一个文件作为代表
            file_stats["key_files"][category] = str(category_files[0].relative_to(Path(data_dir)))
    
    return file_stats


async def map_analysis_to_agent_state(analysis_response: AnalysisResponse, state: AgentState) -> Dict[str, Any]:
    """将AnalysisResponse映射到完整的AgentState"""
    return {
        "messages": state.messages,
        "input": state.input,
        "response": analysis_response.analysis_summary,
        "status": "analysis_complete",
        "analysis_summary": analysis_response.analysis_summary,
        "analysis_insights": analysis_response.analysis_insights,
        "result_files": analysis_response.result_files,
        "quality_metrics": analysis_response.quality_metrics,
        "next_steps": analysis_response.next_steps,
        # 保持执行状态
        "execution_status": state.execution_status,
        "execution_output": state.execution_output,
        "execution_result": state.execution_result,
        "nextflow_config": state.nextflow_config,
    }


async def analysis_node(state: AgentState) -> Dict[str, Any]:
    """
    Analysis节点 - 提取关键指标并让LLM生成智能分析总结
    
    架构：
    1. 从执行输出中提取具体指标
    2. 扫描结果文件统计
    3. 将结构化指标交给LLM分析
    4. LLM生成专业的总结和建议
    """
    
    execution_output = getattr(state, 'execution_output', '')
    execution_status = getattr(state, 'execution_status', '')
    nextflow_config = getattr(state, 'nextflow_config', {})
    
    # 1. 提取各工具的具体指标  
    data_dir = "."  # Docker容器内当前目录就是/data
    extracted_metrics = {
        "nextflow": extract_nextflow_metrics(data_dir),
        "fastp": extract_fastp_metrics(data_dir),
        "star": extract_star_metrics(execution_output), 
        "featurecounts": extract_featurecounts_metrics(data_dir),
        "file_stats": scan_result_files(data_dir)
    }
    
    # 2. 准备配置信息上下文
    analysis_context = {
        "execution_status": execution_status,
        "tools_used": {
            "qc_tool": nextflow_config.get("qc_tool", "unknown"),
            "align_tool": nextflow_config.get("align_tool", "unknown"),
            "quant_tool": nextflow_config.get("quant_tool", "unknown"),
            "genome_version": nextflow_config.get("genome_version", "unknown")
        },
        "extracted_metrics": extracted_metrics
    }
    
    # 3. 让LLM基于具体指标生成专业分析
    llm = get_shared_llm()
    structured_llm = llm.with_structured_output(AnalysisResponse, method="json_mode")
    
    analysis_prompt = f"""你是RNA-seq数据分析专家。请基于以下具体的技术指标生成专业的分析总结报告。

## 分析配置
{json.dumps(analysis_context["tools_used"], ensure_ascii=False, indent=2)}

## 执行状态
{execution_status}

## 具体技术指标
{json.dumps(extracted_metrics, ensure_ascii=False, indent=2)}

请生成专业的RNA-seq分析报告，包含：

1. **analysis_summary**: 基于具体指标的总结(3-4句话，突出关键数值)
   - 如果有比对率，说明比对效果
   - 如果有质控数据，评估数据质量
   - 如果有定量结果，说明基因检出情况

2. **analysis_insights**: 基于数值的专业洞察(每条包含具体数据)
   - 例如："✅ STAR比对成功率达到85.2%，表明样本与参考基因组匹配良好"
   - 例如："📊 featureCounts成功分配了2,345,678条reads到基因特征"

3. **result_files**: 重要结果文件路径
4. **quality_metrics**: 关键质量指标的结构化数据
5. **next_steps**: 基于当前结果的具体建议

要求：
- 使用中文
- 基于实际数值进行评估
- 如果某个指标异常，明确指出问题
- 提供具体可行的改进建议
- 输出JSON格式
"""
    
    try:
        # LLM基于具体指标生成分析
        raw_response = await structured_llm.ainvoke([{"role": "user", "content": analysis_prompt}])
        
        # 确保返回正确的AnalysisResponse类型
        if isinstance(raw_response, dict):
            analysis_response = AnalysisResponse(**raw_response)
        else:
            analysis_response = raw_response
            
        # 确保技术指标被保留
        if not analysis_response.quality_metrics:
            analysis_response.quality_metrics = extracted_metrics
            
    except Exception as e:
        print(f"LLM分析失败: {e}")
        # 备用响应基于提取的指标
        analysis_response = AnalysisResponse(
            analysis_summary=f"分析完成({execution_status})。检测到{len(extracted_metrics)}类技术指标。",
            analysis_insights=[f"工具组合: {analysis_context['tools_used']}"],
            result_files=extracted_metrics.get("file_stats", {}).get("key_files", {}),
            quality_metrics=extracted_metrics,
            next_steps=["检查详细的技术指标进行进一步分析"]
        )
    
    # 4. 保存分析报告到文件
    saved_files = save_analysis_report(analysis_response)
    print(f"✅ 分析报告已保存: {saved_files}")
    
    return await map_analysis_to_agent_state(analysis_response, state)