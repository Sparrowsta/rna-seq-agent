import json
import re
from typing import Dict, Any, List
from pathlib import Path
import pandas as pd

from ..core import get_shared_llm
from ..state import AgentState, AnalysisResponse
from ..prompts import ANALYSIS_NODE_PROMPT
from langgraph.prebuilt import create_react_agent

def create_analysis_agent():
    """创建Analysis节点的智能分析Agent"""
    llm = get_shared_llm()
    
    # 使用create_react_agent但不提供tools，纯推理模式
    agent = create_react_agent(
        model=llm,
        tools=[],  # 空工具列表，纯推理
        prompt=ANALYSIS_NODE_PROMPT,  # 使用集中管理的prompt
        response_format=AnalysisResponse
    )
    return agent

def extract_sample_fastp_metrics(sample_dir: Path) -> Dict[str, Any]:
    """提取单个样本的fastp指标"""
    metrics = {
        'sample_id': sample_dir.name,
        'sequencing_type': 'unknown',
        'total_reads_raw': 0,
        'total_reads_clean': 0,
        'q20_rate': 0.0,
        'q30_rate': 0.0,
        'gc_content': 0.0,
        'duplication_rate': 0.0,
        'adapter_trimmed_reads': 0,
        'reads_filtered_rate': 0.0
    }
    
    json_files = list(sample_dir.glob("*.fastp.json"))
    if not json_files:
        return metrics
    
    try:
        with open(json_files[0], 'r') as f:
            data = json.load(f)
            
        summary = data.get("summary", {})
        filtering = data.get("filtering_result", {})
        duplication = data.get("duplication", {})
        adapter_cutting = data.get("adapter_cutting", {})
        
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})
        
        # 基本信息
        metrics['sample_id'] = sample_dir.name
        sequencing = summary.get("sequencing", "")
        if "paired end" in sequencing.lower():
            metrics['sequencing_type'] = 'paired_end'
        elif "single end" in sequencing.lower():
            metrics['sequencing_type'] = 'single_end'
        
        # 读取数量
        total_raw = before.get("total_reads", 0)
        total_clean = after.get("total_reads", 0)
        metrics['total_reads_raw'] = total_raw
        metrics['total_reads_clean'] = total_clean
        
        # 过滤率
        if total_raw > 0:
            metrics['reads_filtered_rate'] = (total_raw - total_clean) / total_raw
        
        # 质量指标
        metrics['q20_rate'] = after.get("q20_rate", 0.0)
        metrics['q30_rate'] = after.get("q30_rate", 0.0)
        metrics['gc_content'] = after.get("gc_content", 0.0)
        
        # 重复率
        metrics['duplication_rate'] = duplication.get("rate", 0.0)
        
        # 接头修剪
        metrics['adapter_trimmed_reads'] = adapter_cutting.get("adapter_trimmed_reads", 0)
        
    except (json.JSONDecodeError, FileNotFoundError, KeyError) as e:
        print(f"Warning: Failed to parse fastp data for {sample_dir.name}: {e}")
    
    return metrics


def extract_sample_star_metrics(sample_dir: Path) -> Dict[str, Any]:
    """提取单个样本的STAR指标"""
    metrics = {
        'sample_id': sample_dir.name,
        'input_reads': 0,
        'uniquely_mapped_reads': 0,
        'uniquely_mapped_rate': 0.0,
        'multi_mapped_reads': 0,
        'multi_mapped_rate': 0.0,
        'unmapped_reads': 0,
        'unmapped_rate': 0.0,
        'mismatch_rate': 0.0,
        'splice_sites_total': 0
    }
    
    log_files = list(sample_dir.glob("*.star.log"))
    if not log_files:
        return metrics
    
    try:
        with open(log_files[0], 'r', encoding='utf-8') as f:
            log_content = f.read()
        
        # 提取关键统计指标
        patterns = {
            "input_reads": r'Number of input reads \|\s+(\d+)',
            "uniquely_mapped": r'Uniquely mapped reads number \|\s+(\d+)',
            "uniquely_mapped_pct": r'Uniquely mapped reads % \|\s+([\d.]+)%',
            "multi_mapped": r'Number of reads mapped to multiple loci \|\s+(\d+)',
            "multi_mapped_pct": r'% of reads mapped to multiple loci \|\s+([\d.]+)%',
            "unmapped_mismatches": r'Number of reads unmapped: too many mismatches \|\s+(\d+)',
            "unmapped_short": r'Number of reads unmapped: too short \|\s+(\d+)',
            "unmapped_other": r'Number of reads unmapped: other \|\s+(\d+)',
            "mismatch_rate": r'Mismatch rate per base, % \|\s+([\d.]+)%',
            "splice_total": r'Number of splices: Total \|\s+(\d+)'
        }
        
        sample_data = {}
        for key, pattern in patterns.items():
            match = re.search(pattern, log_content)
            if match:
                sample_data[key] = match.group(1)
        
        # 填充指标
        metrics['sample_id'] = sample_dir.name
        
        if "input_reads" in sample_data:
            input_reads = int(sample_data["input_reads"])
            metrics['input_reads'] = input_reads
            
            if "uniquely_mapped" in sample_data:
                uniquely_mapped = int(sample_data["uniquely_mapped"])
                metrics['uniquely_mapped_reads'] = uniquely_mapped
                metrics['uniquely_mapped_rate'] = uniquely_mapped / input_reads if input_reads > 0 else 0.0
            
            if "multi_mapped" in sample_data:
                multi_mapped = int(sample_data["multi_mapped"])
                metrics['multi_mapped_reads'] = multi_mapped
                metrics['multi_mapped_rate'] = multi_mapped / input_reads if input_reads > 0 else 0.0
            
            # 计算总的未比对reads
            unmapped_total = 0
            for key in ["unmapped_mismatches", "unmapped_short", "unmapped_other"]:
                if key in sample_data:
                    unmapped_total += int(sample_data[key])
            
            metrics['unmapped_reads'] = unmapped_total
            metrics['unmapped_rate'] = unmapped_total / input_reads if input_reads > 0 else 0.0
        
        if "mismatch_rate" in sample_data:
            metrics['mismatch_rate'] = float(sample_data["mismatch_rate"]) / 100.0
        
        if "splice_total" in sample_data:
            metrics['splice_sites_total'] = int(sample_data["splice_total"])
            
    except (FileNotFoundError, ValueError, AttributeError) as e:
        print(f"Warning: Failed to parse STAR data for {sample_dir.name}: {e}")
    
    return metrics


def extract_sample_featurecounts_metrics(summary_file: Path) -> Dict[str, Dict[str, Any]]:
    """提取所有样本的featureCounts指标"""
    sample_metrics = {}
    
    if not summary_file.exists():
        return sample_metrics
    
    try:
        # 读取summary文件
        df = pd.read_csv(summary_file, sep='\t')
        
        # 获取样本列（除了Status列）
        sample_columns = [col for col in df.columns if col != 'Status']
        
        for sample_col in sample_columns:
            # 提取样本ID（去掉.bam后缀）
            sample_id = sample_col.replace('.bam', '')
            
            # 计算各项指标
            assigned = df[df['Status'] == 'Assigned'][sample_col].iloc[0] if len(df[df['Status'] == 'Assigned']) > 0 else 0
            total_reads = df[sample_col].sum()
            
            unassigned_unmapped = df[df['Status'] == 'Unassigned_Unmapped'][sample_col].iloc[0] if len(df[df['Status'] == 'Unassigned_Unmapped']) > 0 else 0
            unassigned_multimapping = df[df['Status'] == 'Unassigned_MultiMapping'][sample_col].iloc[0] if len(df[df['Status'] == 'Unassigned_MultiMapping']) > 0 else 0
            unassigned_nofeatures = df[df['Status'] == 'Unassigned_NoFeatures'][sample_col].iloc[0] if len(df[df['Status'] == 'Unassigned_NoFeatures']) > 0 else 0
            unassigned_ambiguity = df[df['Status'] == 'Unassigned_Ambiguity'][sample_col].iloc[0] if len(df[df['Status'] == 'Unassigned_Ambiguity']) > 0 else 0
            
            sample_metrics[sample_id] = {
                'sample_id': sample_id,
                'assigned_reads': int(assigned),
                'assigned_rate': assigned / total_reads if total_reads > 0 else 0.0,
                'unassigned_unmapped': int(unassigned_unmapped),
                'unassigned_multimapping': int(unassigned_multimapping),
                'unassigned_nofeatures': int(unassigned_nofeatures),
                'unassigned_ambiguity': int(unassigned_ambiguity),
                'total_processed_reads': int(total_reads),
                'unassigned_unmapped_rate': unassigned_unmapped / total_reads if total_reads > 0 else 0.0,
                'unassigned_multimapping_rate': unassigned_multimapping / total_reads if total_reads > 0 else 0.0
            }
            
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"Warning: Failed to parse featureCounts data: {e}")
    
    return sample_metrics


def build_sample_quality_dataframe(data_dir: str = ".") -> pd.DataFrame:
    """构建样本质量分析的pandas DataFrame"""
    
    # 收集所有样本的指标
    all_sample_metrics = []
    
    # 1. 提取fastp指标
    fastp_dir = Path(data_dir) / "results" / "fastp"
    fastp_metrics = {}
    if fastp_dir.exists():
        for sample_dir in fastp_dir.iterdir():
            if sample_dir.is_dir():
                metrics = extract_sample_fastp_metrics(sample_dir)
                fastp_metrics[metrics['sample_id']] = metrics
    
    # 2. 提取STAR指标
    star_dir = Path(data_dir) / "results" / "bam"
    star_metrics = {}
    if star_dir.exists():
        for sample_dir in star_dir.iterdir():
            if sample_dir.is_dir():
                metrics = extract_sample_star_metrics(sample_dir)
                star_metrics[metrics['sample_id']] = metrics
    
    # 3. 提取featureCounts指标
    featurecounts_file = Path(data_dir) / "results" / "featurecounts" / "all_samples.counts.txt.summary"
    featurecounts_metrics = extract_sample_featurecounts_metrics(featurecounts_file)
    
    # 4. 合并所有样本的指标
    all_sample_ids = set()
    all_sample_ids.update(fastp_metrics.keys())
    all_sample_ids.update(star_metrics.keys())
    all_sample_ids.update(featurecounts_metrics.keys())
    
    for sample_id in all_sample_ids:
        sample_data = {'sample_id': sample_id}
        
        # 合并fastp指标
        if sample_id in fastp_metrics:
            sample_data.update(fastp_metrics[sample_id])
        
        # 合并STAR指标
        if sample_id in star_metrics:
            star_data = star_metrics[sample_id]
            # 重命名避免冲突
            sample_data.update({
                'star_input_reads': star_data['input_reads'],
                'uniquely_mapped_reads': star_data['uniquely_mapped_reads'],
                'uniquely_mapped_rate': star_data['uniquely_mapped_rate'],
                'multi_mapped_reads': star_data['multi_mapped_reads'],
                'multi_mapped_rate': star_data['multi_mapped_rate'],
                'unmapped_reads': star_data['unmapped_reads'],
                'unmapped_rate': star_data['unmapped_rate'],
                'mismatch_rate': star_data['mismatch_rate'],
                'splice_sites_total': star_data['splice_sites_total']
            })
        
        # 合并featureCounts指标
        if sample_id in featurecounts_metrics:
            fc_data = featurecounts_metrics[sample_id]
            sample_data.update({
                'assigned_reads': fc_data['assigned_reads'],
                'assigned_rate': fc_data['assigned_rate'],
                'unassigned_unmapped': fc_data['unassigned_unmapped'],
                'unassigned_multimapping': fc_data['unassigned_multimapping'],
                'unassigned_nofeatures': fc_data['unassigned_nofeatures'],
                'unassigned_ambiguity': fc_data['unassigned_ambiguity'],
                'fc_total_processed_reads': fc_data['total_processed_reads'],
                'unassigned_unmapped_rate': fc_data['unassigned_unmapped_rate'],
                'unassigned_multimapping_rate': fc_data['unassigned_multimapping_rate']
            })
        
        all_sample_metrics.append(sample_data)
    
    # 5. 创建DataFrame
    if all_sample_metrics:
        df = pd.DataFrame(all_sample_metrics)
        
        # 添加质量评估
        df = add_quality_assessment(df)
        
        return df
    else:
        # 返回空DataFrame但包含预期的列
        return pd.DataFrame(columns=[
            'sample_id', 'sequencing_type', 'total_reads_raw', 'total_reads_clean',
            'q20_rate', 'q30_rate', 'gc_content', 'duplication_rate',
            'uniquely_mapped_rate', 'multi_mapped_rate', 'unmapped_rate',
            'assigned_rate', 'quality_score', 'quality_status'
        ])


def add_quality_assessment(df: pd.DataFrame) -> pd.DataFrame:
    """为DataFrame添加质量评估"""
    
    # 质量阈值
    thresholds = {
        'min_mapping_rate': 0.70,
        'max_duplication_rate': 0.30,
        'min_assignment_rate': 0.60,
        'min_q30_rate': 0.80,
        'gc_content_range': (0.40, 0.60),
        'max_mismatch_rate': 0.05
    }
    
    # 计算质量分数和状态
    quality_scores = []
    quality_statuses = []
    quality_flags = []
    
    for _, row in df.iterrows():
        score = 100  # 起始分数
        flags = []
        
        # 检查比对率
        if 'uniquely_mapped_rate' in row and pd.notna(row['uniquely_mapped_rate']):
            if row['uniquely_mapped_rate'] < thresholds['min_mapping_rate']:
                score -= 20
                flags.append('low_mapping_rate')
        
        # 检查重复率
        if 'duplication_rate' in row and pd.notna(row['duplication_rate']):
            if row['duplication_rate'] > thresholds['max_duplication_rate']:
                score -= 15
                flags.append('high_duplication')
        
        # 检查基因分配率
        if 'assigned_rate' in row and pd.notna(row['assigned_rate']):
            if row['assigned_rate'] < thresholds['min_assignment_rate']:
                score -= 20
                flags.append('low_assignment_rate')
        
        # 检查Q30率
        if 'q30_rate' in row and pd.notna(row['q30_rate']):
            if row['q30_rate'] < thresholds['min_q30_rate']:
                score -= 15
                flags.append('low_q30_rate')
        
        # 检查GC含量
        if 'gc_content' in row and pd.notna(row['gc_content']):
            gc_min, gc_max = thresholds['gc_content_range']
            if row['gc_content'] < gc_min or row['gc_content'] > gc_max:
                score -= 10
                flags.append('abnormal_gc_content')
        
        # 检查错配率
        if 'mismatch_rate' in row and pd.notna(row['mismatch_rate']):
            if row['mismatch_rate'] > thresholds['max_mismatch_rate']:
                score -= 10
                flags.append('high_mismatch_rate')
        
        # 确保分数不低于0
        score = max(0, score)
        
        # 确定质量状态
        if score >= 80:
            status = 'PASS'
        elif score >= 60:
            status = 'WARNING'
        else:
            status = 'FAIL'
        
        quality_scores.append(score)
        quality_statuses.append(status)
        quality_flags.append(flags)
    
    # 添加质量评估列
    df['quality_score'] = quality_scores
    df['quality_status'] = quality_statuses
    df['quality_flags'] = quality_flags
    
    return df








def save_analysis_report(analysis_response: AnalysisResponse, data_dir: str = ".") -> Dict[str, str]:
    """保存分析报告到基于时间戳的归档文件夹"""
    import datetime
    import shutil
    
    reports_path = Path(data_dir) / "reports"
    reports_path.mkdir(exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 创建时间戳文件夹
    archive_folder = reports_path / timestamp
    archive_folder.mkdir(exist_ok=True)
    
    saved_files = {}
    
    try:
        # 1. 保存JSON格式的完整数据
        json_file = archive_folder / "analysis_report.json"
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
        md_file = archive_folder / "analysis_summary.md"
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
        
        # 3. 复制runtime_config.json到归档文件夹
        runtime_config_source = Path(data_dir) / "config" / "runtime_config.json"
        if runtime_config_source.exists():
            runtime_config_dest = archive_folder / "runtime_config.json"
            shutil.copy2(runtime_config_source, runtime_config_dest)
            saved_files["runtime_config"] = str(runtime_config_dest.relative_to(Path(data_dir)))
        
        # 4. 创建指向最新归档文件夹的符号链接
        latest_link = reports_path / "latest"
        if latest_link.is_symlink() or latest_link.exists():
            latest_link.unlink()
        latest_link.symlink_to(archive_folder.name)
        saved_files["latest_folder"] = str(latest_link.relative_to(Path(data_dir)))
        
        # 5. 创建执行日志文件占位符(供后续使用)
        log_file = archive_folder / "execution_log.txt"
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write(f"# RNA-seq 分析执行日志\n")
            f.write(f"生成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"分析ID: {timestamp}\n\n")
            f.write("详细的执行日志将在后续版本中添加。\n")
        saved_files["execution_log"] = str(log_file.relative_to(Path(data_dir)))
        
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
        "status": "user_communication",
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
    Analysis节点 - 基于样本级别的表格化数据分析
    
    新架构：
    1. 构建样本质量分析DataFrame
    2. 进行质量异常检测和评估
    3. 生成表格化的质量报告
    4. LLM基于表格数据生成专业分析
    """
    
    execution_output = getattr(state, 'execution_output', '')
    execution_status = getattr(state, 'execution_status', '')
    nextflow_config = getattr(state, 'nextflow_config', {})
    
    # 1. 构建样本质量分析DataFrame
    data_dir = "."  # Docker容器内当前目录就是/data
    
    try:
        quality_df = build_sample_quality_dataframe(data_dir)
        print(f"✅ 构建了包含 {len(quality_df)} 个样本的质量分析表格")
        
        # 2. 生成完整的聚合分析
        comprehensive_summary = generate_comprehensive_summary(quality_df)
        
        # 3. 保留原有的文件统计（用于兼容性）
        file_stats = scan_result_files(data_dir)
        
        # 4. 构建新的分析上下文
        analysis_context = {
            "tools_used": {
                "qc_tool": nextflow_config.get("qc_tool", "fastp"),
                "align_tool": nextflow_config.get("align_tool", "STAR"),
                "quant_tool": nextflow_config.get("quant_tool", "featureCounts"),
                "genome_version": nextflow_config.get("genome_version", "unknown")
            },
            "sample_count": len(quality_df),
            "comprehensive_summary": comprehensive_summary,
            "sample_quality_table": quality_df.to_dict('records') if not quality_df.empty else [],
            "file_stats": file_stats
        }
        
    except Exception as e:
        print(f"Warning: 表格化分析失败，回退到空DataFrame: {e}")
        # 回退到空DataFrame，保持一致的数据结构
        quality_df = pd.DataFrame()
        comprehensive_summary = {"total_samples": 0}
        file_stats = scan_result_files(data_dir)
        
        analysis_context = {
            "tools_used": {
                "qc_tool": nextflow_config.get("qc_tool", "unknown"),
                "align_tool": nextflow_config.get("align_tool", "unknown"),
                "quant_tool": nextflow_config.get("quant_tool", "unknown"),
                "genome_version": nextflow_config.get("genome_version", "unknown")
            },
            "sample_count": 0,
            "comprehensive_summary": comprehensive_summary,
            "sample_quality_table": [],
            "file_stats": file_stats
        }
    
    # 5. 让LLM基于表格化数据生成专业分析
    agent_executor = create_analysis_agent()
    
    analysis_prompt = f"""你是RNA-seq数据分析专家。请基于以下样本级别的质量分析数据生成专业的分析总结报告。

## 分析配置
{json.dumps(analysis_context["tools_used"], ensure_ascii=False, indent=2)}

## 样本质量分析数据
{json.dumps(analysis_context, ensure_ascii=False, indent=2)}

请生成专业的RNA-seq分析报告，包含：

1. **analysis_summary**: 基于样本级别数据的总结(3-4句话，突出关键发现)
   - 总体样本数量和质量状态分布
   - 关键质量指标的表现（比对率、分配率等）
   - 是否发现异常样本及其特征

2. **analysis_insights**: 基于样本数据的专业洞察(每条包含具体数据)
   - 例如："✅ 3个样本中有2个达到PASS标准，平均比对率为85.2%"
   - 例如："⚠️ 样本SRR123456的比对率仅为15.3%，可能存在样本质量问题"
   - 例如："📊 所有样本的基因分配率均超过60%，定量结果可靠"

3. **result_files**: 重要结果文件路径
4. **quality_metrics**: 样本质量分析的结构化数据
5. **next_steps**: 基于样本质量评估的具体建议

要求：
- 使用中文
- 重点关注样本级别的质量差异
- 明确指出质量异常的样本
- 提供针对性的改进建议
- 输出JSON格式
"""
    
    try:
        # LLM基于表格化数据生成分析
        messages_input = {"messages": [{"role": "user", "content": analysis_prompt}]}
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")
        
        if structured_response:
            analysis_response = structured_response
        else:
            raise Exception("Agent未返回预期的结构化响应")
            
        # 确保质量指标被保留
        if not analysis_response.quality_metrics:
            analysis_response.quality_metrics = analysis_context
            
    except Exception as e:
        print(f"LLM分析失败: {e}")
        # 备用响应基于表格化数据
        if 'comprehensive_summary' in analysis_context:
            summary_text = f"完成了{analysis_context['sample_count']}个样本的质量分析。"
            if analysis_context['comprehensive_summary'].get('fail_samples', 0) > 0:
                summary_text += f"发现{analysis_context['comprehensive_summary']['fail_samples']}个质量异常样本。"
        else:
            summary_text = "分析完成，使用传统聚合方法。"
            
        analysis_response = AnalysisResponse(
            analysis_summary=summary_text,
            analysis_insights=["已完成样本级别的质量分析"],
            result_files=analysis_context.get("file_stats", {}).get("key_files", {}),
            quality_metrics=analysis_context,
            next_steps=["检查样本质量表格进行详细分析"]
        )
    
    # 6. 保存分析报告到文件
    saved_files = save_analysis_report(analysis_response)
    print(f"✅ 分析报告已保存: {saved_files}")
    
    return await map_analysis_to_agent_state(analysis_response, state)


def generate_comprehensive_summary(df: pd.DataFrame) -> Dict[str, Any]:
    """从DataFrame生成完整的聚合分析 - 包括质量评估和工具统计"""
    if df.empty:
        return {"total_samples": 0}
    
    summary: Dict[str, Any] = {
        "total_samples": len(df),
    }
    
    # ========== 质量评估聚合 ==========
    if 'quality_status' in df.columns:
        summary.update({
            "pass_samples": len(df[df['quality_status'] == 'PASS']),
            "warning_samples": len(df[df['quality_status'] == 'WARNING']),
            "fail_samples": len(df[df['quality_status'] == 'FAIL']),
        })
        
        # 异常样本识别
        fail_samples = df[df['quality_status'] == 'FAIL']['sample_id'].tolist()
        warning_samples = df[df['quality_status'] == 'WARNING']['sample_id'].tolist()
        summary["fail_sample_ids"] = fail_samples
        summary["warning_sample_ids"] = warning_samples
    
    # 关键指标统计（均值、最值）
    numeric_columns = ['uniquely_mapped_rate', 'assigned_rate', 'q30_rate', 'duplication_rate']
    for col in numeric_columns:
        if col in df.columns and df[col].notna().any():
            summary[f"{col}_mean"] = float(df[col].mean())
            summary[f"{col}_min"] = float(df[col].min())
            summary[f"{col}_max"] = float(df[col].max())
    
    # ========== 工具数据聚合 ==========
    # fastp聚合指标
    if 'total_reads_raw' in df.columns and 'total_reads_clean' in df.columns:
        total_raw = df['total_reads_raw'].sum()
        total_clean = df['total_reads_clean'].sum()
        
        summary["fastp_summary"] = {
            "total_before_reads": int(total_raw),
            "total_after_reads": int(total_clean),
            "filtering_rate": f"{((total_raw - total_clean) / total_raw * 100):.1f}%" if total_raw > 0 else "0%",
            "average_q30_rate": f"{df['q30_rate'].mean() * 100:.1f}%" if 'q30_rate' in df.columns else "0%",
            "average_duplication_rate": f"{df['duplication_rate'].mean() * 100:.1f}%" if 'duplication_rate' in df.columns else "0%"
        }
    
    # STAR聚合指标
    if 'star_input_reads' in df.columns and 'uniquely_mapped_reads' in df.columns:
        total_input = df['star_input_reads'].sum()
        total_uniquely = df['uniquely_mapped_reads'].sum()
        
        summary["star_summary"] = {
            "total_input_reads": int(total_input),
            "total_uniquely_mapped": int(total_uniquely),
            "average_mapping_rate": f"{df['uniquely_mapped_rate'].mean() * 100:.1f}%" if 'uniquely_mapped_rate' in df.columns else "0%"
        }
    
    # featureCounts聚合指标
    if 'assigned_reads' in df.columns and 'fc_total_processed_reads' in df.columns:
        total_assigned = df['assigned_reads'].sum()
        total_processed = df['fc_total_processed_reads'].sum()
        
        summary["featurecounts_summary"] = {
            "total_assigned_reads": int(total_assigned),
            "total_processed_reads": int(total_processed),
            "average_assignment_rate": f"{df['assigned_rate'].mean() * 100:.1f}%" if 'assigned_rate' in df.columns else "0%"
        }
    
    return summary
