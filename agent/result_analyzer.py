"""
结果分析和总结模块
遵循单一职责原则：专门处理RNA-seq分析结果的收集、解析和总结
"""

import os
import json
import logging
import re
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass, asdict
from datetime import datetime

# 配置日志
logger = logging.getLogger(__name__)

# ============================================================================
# 结果数据模型 - 遵循数据类模式
# ============================================================================

@dataclass
class FastpResult:
    """FastP质量控制结果"""
    sample_id: str
    total_reads_before: int = 0
    total_reads_after: int = 0
    total_bases_before: int = 0
    total_bases_after: int = 0
    q20_rate_before: float = 0.0
    q20_rate_after: float = 0.0
    q30_rate_before: float = 0.0
    q30_rate_after: float = 0.0
    gc_content_before: float = 0.0
    gc_content_after: float = 0.0
    duplication_rate: float = 0.0
    
    def get_summary(self) -> Dict[str, Any]:
        """获取质控结果摘要"""
        return {
            "样本ID": self.sample_id,
            "质控前读数": f"{self.total_reads_before:,}",
            "质控后读数": f"{self.total_reads_after:,}",
            "读数保留率": f"{(self.total_reads_after/self.total_reads_before*100):.1f}%" if self.total_reads_before > 0 else "N/A",
            "Q30质量": f"{self.q30_rate_after:.1f}%",
            "GC含量": f"{self.gc_content_after:.1f}%",
            "重复率": f"{self.duplication_rate:.1f}%"
        }

@dataclass
class StarResult:
    """STAR比对结果"""
    sample_id: str
    input_reads: int = 0
    uniquely_mapped: int = 0
    uniquely_mapped_percent: float = 0.0
    multi_mapped: int = 0
    multi_mapped_percent: float = 0.0
    unmapped_too_many_mismatches: int = 0
    unmapped_too_short: int = 0
    unmapped_other: int = 0
    
    def get_summary(self) -> Dict[str, Any]:
        """获取比对结果摘要"""
        total_mapped = self.uniquely_mapped + self.multi_mapped
        total_mapped_percent = (total_mapped / self.input_reads * 100) if self.input_reads > 0 else 0
        
        return {
            "样本ID": self.sample_id,
            "输入读数": f"{self.input_reads:,}",
            "唯一比对": f"{self.uniquely_mapped:,} ({self.uniquely_mapped_percent:.1f}%)",
            "多重比对": f"{self.multi_mapped:,} ({self.multi_mapped_percent:.1f}%)",
            "总比对率": f"{total_mapped_percent:.1f}%",
            "未比对读数": f"{self.input_reads - total_mapped:,}"
        }

@dataclass
class FeatureCountsResult:
    """FeatureCounts定量结果"""
    total_reads: int = 0
    assigned_reads: int = 0
    unassigned_unmapped: int = 0
    unassigned_read_type: int = 0
    unassigned_singleton: int = 0
    unassigned_mapping_quality: int = 0
    unassigned_chimera: int = 0
    unassigned_secondary: int = 0
    unassigned_nonjunction: int = 0
    unassigned_duplicate: int = 0
    gene_count: int = 0
    
    def get_summary(self) -> Dict[str, Any]:
        """获取定量结果摘要"""
        assignment_rate = (self.assigned_reads / self.total_reads * 100) if self.total_reads > 0 else 0
        
        return {
            "总读数": f"{self.total_reads:,}",
            "成功分配": f"{self.assigned_reads:,} ({assignment_rate:.1f}%)",
            "检测基因数": f"{self.gene_count:,}",
            "未分配读数": f"{self.total_reads - self.assigned_reads:,}",
            "分配效率": f"{assignment_rate:.1f}%"
        }

@dataclass
class AnalysisResult:
    """完整分析结果"""
    analysis_id: str
    start_time: str
    end_time: str
    status: str
    fastp_results: List[FastpResult]
    star_results: List[StarResult]
    featurecounts_result: Optional[FeatureCountsResult]
    output_files: Dict[str, List[str]]
    
    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式"""
        return asdict(self)

# ============================================================================
# 结果解析器 - 遵循策略模式
# ============================================================================

class FastpResultParser:
    """
    FastP结果解析器
    
    遵循单一职责原则：专门解析FastP的JSON输出
    """
    
    @staticmethod
    def parse_json_file(json_file: str) -> Optional[FastpResult]:
        """解析FastP的JSON结果文件"""
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # 提取样本ID
            sample_id = Path(json_file).stem.replace('.fastp', '')
            
            # 解析统计数据
            summary = data.get('summary', {})
            before_filtering = summary.get('before_filtering', {})
            after_filtering = summary.get('after_filtering', {})
            
            result = FastpResult(
                sample_id=sample_id,
                total_reads_before=before_filtering.get('total_reads', 0),
                total_reads_after=after_filtering.get('total_reads', 0),
                total_bases_before=before_filtering.get('total_bases', 0),
                total_bases_after=after_filtering.get('total_bases', 0),
                q20_rate_before=before_filtering.get('q20_rate', 0.0) * 100,
                q20_rate_after=after_filtering.get('q20_rate', 0.0) * 100,
                q30_rate_before=before_filtering.get('q30_rate', 0.0) * 100,
                q30_rate_after=after_filtering.get('q30_rate', 0.0) * 100,
                gc_content_before=before_filtering.get('gc_content', 0.0) * 100,
                gc_content_after=after_filtering.get('gc_content', 0.0) * 100,
                duplication_rate=data.get('duplication', {}).get('rate', 0.0) * 100
            )
            
            logger.info(f"成功解析FastP结果: {sample_id}")
            return result
        
        except Exception as e:
            logger.error(f"解析FastP结果文件失败 {json_file}: {e}")
            return None

class StarResultParser:
    """
    STAR结果解析器
    
    遵循单一职责原则：专门解析STAR的Log.final.out文件
    """
    
    @staticmethod
    def parse_log_file(log_file: str) -> Optional[StarResult]:
        """解析STAR的Log.final.out文件"""
        try:
            with open(log_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 提取样本ID
            sample_id = Path(log_file).parent.name
            
            # 使用正则表达式提取统计信息
            patterns = {
                'input_reads': r'Number of input reads \|\s+(\d+)',
                'uniquely_mapped': r'Uniquely mapped reads number \|\s+(\d+)',
                'uniquely_mapped_percent': r'Uniquely mapped reads % \|\s+([\d.]+)%',
                'multi_mapped': r'Number of reads mapped to multiple loci \|\s+(\d+)',
                'multi_mapped_percent': r'% of reads mapped to multiple loci \|\s+([\d.]+)%',
                'unmapped_too_many_mismatches': r'% of reads unmapped: too many mismatches \|\s+([\d.]+)%',
                'unmapped_too_short': r'% of reads unmapped: too short \|\s+([\d.]+)%',
                'unmapped_other': r'% of reads unmapped: other \|\s+([\d.]+)%'
            }
            
            extracted_data = {}
            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    value = match.group(1)
                    if 'percent' in key:
                        extracted_data[key] = float(value)
                    else:
                        extracted_data[key] = int(value)
                else:
                    extracted_data[key] = 0
            
            result = StarResult(
                sample_id=sample_id,
                input_reads=extracted_data.get('input_reads', 0),
                uniquely_mapped=extracted_data.get('uniquely_mapped', 0),
                uniquely_mapped_percent=extracted_data.get('uniquely_mapped_percent', 0.0),
                multi_mapped=extracted_data.get('multi_mapped', 0),
                multi_mapped_percent=extracted_data.get('multi_mapped_percent', 0.0),
                unmapped_too_many_mismatches=extracted_data.get('unmapped_too_many_mismatches', 0),
                unmapped_too_short=extracted_data.get('unmapped_too_short', 0),
                unmapped_other=extracted_data.get('unmapped_other', 0)
            )
            
            logger.info(f"成功解析STAR结果: {sample_id}")
            return result
        
        except Exception as e:
            logger.error(f"解析STAR结果文件失败 {log_file}: {e}")
            return None

class FeatureCountsResultParser:
    """
    FeatureCounts结果解析器
    
    遵循单一职责原则：专门解析FeatureCounts的summary文件
    """
    
    @staticmethod
    def parse_summary_file(summary_file: str) -> Optional[FeatureCountsResult]:
        """解析FeatureCounts的summary文件"""
        try:
            with open(summary_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            # 解析统计信息
            stats = {}
            for line in lines[1:]:  # 跳过标题行
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    category = parts[0]
                    count = int(parts[1])
                    stats[category] = count
            
            # 计算基因数量（从counts文件）
            gene_count = 0
            counts_file = summary_file.replace('.summary', '')
            if os.path.exists(counts_file):
                with open(counts_file, 'r') as f:
                    gene_count = sum(1 for line in f) - 1  # 减去标题行
            
            result = FeatureCountsResult(
                total_reads=sum(stats.values()),
                assigned_reads=stats.get('Assigned', 0),
                unassigned_unmapped=stats.get('Unassigned_Unmapped', 0),
                unassigned_read_type=stats.get('Unassigned_Read_Type', 0),
                unassigned_singleton=stats.get('Unassigned_Singleton', 0),
                unassigned_mapping_quality=stats.get('Unassigned_MappingQuality', 0),
                unassigned_chimera=stats.get('Unassigned_Chimera', 0),
                unassigned_secondary=stats.get('Unassigned_Secondary', 0),
                unassigned_nonjunction=stats.get('Unassigned_NoFeatures', 0),
                unassigned_duplicate=stats.get('Unassigned_Duplicate', 0),
                gene_count=gene_count
            )
            
            logger.info("成功解析FeatureCounts结果")
            return result
        
        except Exception as e:
            logger.error(f"解析FeatureCounts结果文件失败 {summary_file}: {e}")
            return None

# ============================================================================
# 结果收集器 - 遵循收集器模式
# ============================================================================

class ResultCollector:
    """
    结果收集器
    
    遵循单一职责原则：专门收集和整合各种分析结果
    """
    
    def __init__(self, results_dir: str):
        self.results_dir = Path(results_dir)
        self.fastp_parser = FastpResultParser()
        self.star_parser = StarResultParser()
        self.featurecounts_parser = FeatureCountsResultParser()
    
    def collect_fastp_results(self) -> List[FastpResult]:
        """收集所有FastP结果"""
        results = []
        fastp_dir = self.results_dir / "fastp"
        
        if not fastp_dir.exists():
            logger.warning(f"FastP结果目录不存在: {fastp_dir}")
            return results
        
        # 查找所有FastP JSON文件
        json_files = list(fastp_dir.rglob("*.fastp.json"))
        
        for json_file in json_files:
            result = self.fastp_parser.parse_json_file(str(json_file))
            if result:
                results.append(result)
        
        logger.info(f"收集到 {len(results)} 个FastP结果")
        return results
    
    def collect_star_results(self) -> List[StarResult]:
        """收集所有STAR结果"""
        results = []
        star_dir = self.results_dir / "bam"
        
        if not star_dir.exists():
            logger.warning(f"STAR结果目录不存在: {star_dir}")
            return results
        
        # 查找所有STAR日志文件
        log_files = list(star_dir.rglob("Log.final.out"))
        
        for log_file in log_files:
            result = self.star_parser.parse_log_file(str(log_file))
            if result:
                results.append(result)
        
        logger.info(f"收集到 {len(results)} 个STAR结果")
        return results
    
    def collect_featurecounts_results(self) -> Optional[FeatureCountsResult]:
        """收集FeatureCounts结果"""
        featurecounts_dir = self.results_dir / "featurecounts"
        
        if not featurecounts_dir.exists():
            logger.warning(f"FeatureCounts结果目录不存在: {featurecounts_dir}")
            return None
        
        # 查找summary文件
        summary_files = list(featurecounts_dir.glob("*.summary"))
        
        if not summary_files:
            logger.warning("未找到FeatureCounts summary文件")
            return None
        
        # 使用第一个找到的summary文件
        summary_file = summary_files[0]
        result = self.featurecounts_parser.parse_summary_file(str(summary_file))
        
        if result:
            logger.info("成功收集FeatureCounts结果")
        
        return result
    
    def collect_output_files(self) -> Dict[str, List[str]]:
        """收集所有输出文件"""
        output_files = {
            "fastp": [],
            "star": [],
            "featurecounts": [],
            "logs": []
        }
        
        try:
            # FastP文件
            fastp_dir = self.results_dir / "fastp"
            if fastp_dir.exists():
                output_files["fastp"] = [str(f) for f in fastp_dir.rglob("*") if f.is_file()]
            
            # STAR文件
            star_dir = self.results_dir / "bam"
            if star_dir.exists():
                output_files["star"] = [str(f) for f in star_dir.rglob("*") if f.is_file()]
            
            # FeatureCounts文件
            fc_dir = self.results_dir / "featurecounts"
            if fc_dir.exists():
                output_files["featurecounts"] = [str(f) for f in fc_dir.rglob("*") if f.is_file()]
            
            # 日志文件
            logs_dir = self.results_dir.parent / "logs"
            if logs_dir.exists():
                output_files["logs"] = [str(f) for f in logs_dir.rglob("*") if f.is_file()]
        
        except Exception as e:
            logger.error(f"收集输出文件时出错: {e}")
        
        return output_files
    
    def collect_all_results(self, analysis_id: str = None) -> AnalysisResult:
        """收集所有分析结果"""
        if analysis_id is None:
            analysis_id = f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        logger.info(f"开始收集分析结果: {analysis_id}")
        
        # 收集各类结果
        fastp_results = self.collect_fastp_results()
        star_results = self.collect_star_results()
        featurecounts_result = self.collect_featurecounts_results()
        output_files = self.collect_output_files()
        
        # 创建完整结果对象
        result = AnalysisResult(
            analysis_id=analysis_id,
            start_time=datetime.now().isoformat(),
            end_time=datetime.now().isoformat(),
            status="completed",
            fastp_results=fastp_results,
            star_results=star_results,
            featurecounts_result=featurecounts_result,
            output_files=output_files
        )
        
        logger.info(f"结果收集完成: {len(fastp_results)} FastP, {len(star_results)} STAR, {'1' if featurecounts_result else '0'} FeatureCounts")
        
        return result

# ============================================================================
# 结果报告生成器 - 遵循模板方法模式
# ============================================================================

class ResultReportGenerator:
    """
    结果报告生成器
    
    遵循模板方法模式：提供标准的报告生成流程
    """
    
    def __init__(self):
        self.report_template = self._load_report_template()
    
    def _load_report_template(self) -> str:
        """加载报告模板"""
        return """
# RNA-seq分析结果报告

## 分析概览
- **分析ID**: {analysis_id}
- **开始时间**: {start_time}
- **结束时间**: {end_time}
- **分析状态**: {status}

## 质量控制结果 (FastP)
{fastp_summary}

## 序列比对结果 (STAR)
{star_summary}

## 基因定量结果 (FeatureCounts)
{featurecounts_summary}

## 输出文件统计
{output_files_summary}

## 分析建议
{recommendations}

---
*报告生成时间: {report_time}*
        """.strip()
    
    def generate_fastp_summary(self, fastp_results: List[FastpResult]) -> str:
        """生成FastP结果摘要"""
        if not fastp_results:
            return "未找到FastP质量控制结果。"
        
        summary_lines = ["| 样本ID | 质控前读数 | 质控后读数 | 保留率 | Q30质量 | GC含量 |"]
        summary_lines.append("|--------|------------|------------|--------|---------|--------|")
        
        for result in fastp_results:
            summary = result.get_summary()
            line = f"| {summary['样本ID']} | {summary['质控前读数']} | {summary['质控后读数']} | {summary['读数保留率']} | {summary['Q30质量']} | {summary['GC含量']} |"
            summary_lines.append(line)
        
        return "\n".join(summary_lines)
    
    def generate_star_summary(self, star_results: List[StarResult]) -> str:
        """生成STAR结果摘要"""
        if not star_results:
            return "未找到STAR比对结果。"
        
        summary_lines = ["| 样本ID | 输入读数 | 唯一比对 | 多重比对 | 总比对率 |"]
        summary_lines.append("|--------|----------|----------|----------|----------|")
        
        for result in star_results:
            summary = result.get_summary()
            line = f"| {summary['样本ID']} | {summary['输入读数']} | {summary['唯一比对']} | {summary['多重比对']} | {summary['总比对率']} |"
            summary_lines.append(line)
        
        return "\n".join(summary_lines)
    
    def generate_featurecounts_summary(self, fc_result: Optional[FeatureCountsResult]) -> str:
        """生成FeatureCounts结果摘要"""
        if not fc_result:
            return "未找到FeatureCounts定量结果。"
        
        summary = fc_result.get_summary()
        return f"""
- **总读数**: {summary['总读数']}
- **成功分配**: {summary['成功分配']}
- **检测基因数**: {summary['检测基因数']}
- **分配效率**: {summary['分配效率']}
        """.strip()
    
    def generate_output_files_summary(self, output_files: Dict[str, List[str]]) -> str:
        """生成输出文件摘要"""
        summary_lines = []
        
        for category, files in output_files.items():
            if files:
                summary_lines.append(f"- **{category.upper()}**: {len(files)} 个文件")
        
        return "\n".join(summary_lines) if summary_lines else "未找到输出文件。"
    
    def generate_recommendations(self, result: AnalysisResult) -> str:
        """生成分析建议"""
        recommendations = []
        
        # 基于FastP结果的建议
        if result.fastp_results:
            avg_q30 = sum(r.q30_rate_after for r in result.fastp_results) / len(result.fastp_results)
            if avg_q30 < 80:
                recommendations.append("- 数据质量较低，建议检查测序质量或调整质控参数")
        
        # 基于STAR结果的建议
        if result.star_results:
            avg_mapping = sum(r.uniquely_mapped_percent for r in result.star_results) / len(result.star_results)
            if avg_mapping < 70:
                recommendations.append("- 比对率较低，建议检查参考基因组或测序数据质量")
        
        # 基于FeatureCounts结果的建议
        if result.featurecounts_result:
            assignment_rate = (result.featurecounts_result.assigned_reads / 
                             result.featurecounts_result.total_reads * 100) if result.featurecounts_result.total_reads > 0 else 0
            if assignment_rate < 60:
                recommendations.append("- 基因分配率较低，建议检查GTF注释文件或比对参数")
        
        return "\n".join(recommendations) if recommendations else "分析结果正常，无特殊建议。"
    
    def generate_report(self, result: AnalysisResult) -> str:
        """生成完整报告"""
        try:
            report = self.report_template.format(
                analysis_id=result.analysis_id,
                start_time=result.start_time,
                end_time=result.end_time,
                status=result.status,
                fastp_summary=self.generate_fastp_summary(result.fastp_results),
                star_summary=self.generate_star_summary(result.star_results),
                featurecounts_summary=self.generate_featurecounts_summary(result.featurecounts_result),
                output_files_summary=self.generate_output_files_summary(result.output_files),
                recommendations=self.generate_recommendations(result),
                report_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            )
            
            return report
        
        except Exception as e:
            logger.error(f"生成报告时出错: {e}")
            return f"报告生成失败: {str(e)}"
    
    def save_report(self, result: AnalysisResult, output_path: str) -> bool:
        """保存报告到文件"""
        try:
            report = self.generate_report(result)
            
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(report)
            
            logger.info(f"报告已保存到: {output_path}")
            return True
        
        except Exception as e:
            logger.error(f"保存报告失败: {e}")
            return False

# ============================================================================
# 便捷函数 - 遵循DRY原则
# ============================================================================

def analyze_results(results_dir: str, analysis_id: str = None) -> AnalysisResult:
    """分析结果目录中的所有结果"""
    collector = ResultCollector(results_dir)
    return collector.collect_all_results(analysis_id)

def generate_analysis_report(results_dir: str, output_path: str = None, analysis_id: str = None) -> str:
    """生成分析报告"""
    # 收集结果
    result = analyze_results(results_dir, analysis_id)
    
    # 生成报告
    generator = ResultReportGenerator()
    report = generator.generate_report(result)
    
    # 保存报告（如果指定了输出路径）
    if output_path:
        generator.save_report(result, output_path)
    
    return report

def get_analysis_summary(results_dir: str) -> Dict[str, Any]:
    """获取分析结果摘要"""
    result = analyze_results(results_dir)
    
    return {
        "analysis_id": result.analysis_id,
        "status": result.status,
        "sample_count": len(result.fastp_results),
        "fastp_completed": len(result.fastp_results) > 0,
        "star_completed": len(result.star_results) > 0,
        "featurecounts_completed": result.featurecounts_result is not None,
        "total_output_files": sum(len(files) for files in result.output_files.values())
    }