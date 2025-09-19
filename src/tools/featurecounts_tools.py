"""
RNA-seq智能分析助手 - FeatureCounts定量工具

包含：
- run_nextflow_featurecounts: 执行FeatureCounts定量流程
- parse_featurecounts_metrics: 解析FeatureCounts结果
"""

import json
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List

import pandas as pd

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config, Settings
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.featurecounts")


@tool
def run_nextflow_featurecounts(
    featurecounts_params: Dict[str, Any],
    star_results: Dict[str, Any],
    gtf_path: str,
    results_timestamp: Optional[str] = None,
    hisat2_results: Optional[Dict[str, Any]] = None,
    resource_config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """执行Nextflow FeatureCounts定量流程

    Args:
        featurecounts_params: FeatureCounts参数字典（可包含 -T/-s/-p/-M 等风格键）
        star_results: STAR节点结果（可为空），需包含 per_sample_outputs 中的 BAM 路径
        gtf_path: GTF注释文件的绝对路径（从节点预处理提供）
        results_timestamp: 可选的时间戳，优先用于结果目录
        hisat2_results: HISAT2节点结果（可为空），与 star_results 二选一

    Returns:
        执行结果字典，包含状态、输出路径、样本输出等
    """
    try:
        tools_config = get_tools_config()

        # 校验依赖输入
        # 选择可用的比对结果（STAR/HISAT2）
        align_results = None
        if star_results and star_results.get("success"):
            align_results = star_results
        elif hisat2_results and hisat2_results.get("success"):
            align_results = hisat2_results
        else:
            return {"success": False, "error": "比对结果无效（STAR/HISAT2），无法执行FeatureCounts"}

        per_sample = align_results.get("per_sample_outputs") or []
        if not per_sample:
            return {"success": False, "error": "比对结果缺少 per_sample_outputs 信息"}


        # 直接使用传入的GTF路径
        gtf_file = gtf_path

        if not gtf_file:
            return {
                "success": False,
                "error": "GTF路径为空，无法执行FeatureCounts"
            }

        # 验证GTF文件是否存在
        if not Path(gtf_file).exists():
            return {
                "success": False,
                "error": f"GTF文件不存在: {gtf_file}",
            }

        # 运行根目录（results_dir）：复用比对步骤的 results_dir，保持同一运行根目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        run_root = Path(align_results.get("results_dir") or tools_config.results_dir / f"{timestamp}")
        results_dir = run_root
        # 统一 Nextflow 工作目录到 /data/work，使用运行ID区分
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"featurecounts_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        (results_dir / "featurecounts").mkdir(parents=True, exist_ok=True)
        logger.info(f"FeatureCounts启动：bam={len(per_sample)} results={results_dir} work={work_dir}")

        # 构建 Nextflow 参数（与 featurecounts.nf 对齐）
        # 将对齐结果（STAR/HISAT2）产出的 BAM 列表转换为 JSON 字符串（Nextflow 端会 echo 后再解析）
        bam_entries = []
        for item in per_sample:
            sample_id = item.get("sample_id") or "sample"
            bam = item.get("aligned_bam") or item.get("bam")
            if not bam:
                continue
            bam_norm = bam
            if not Path(bam_norm).exists():
                return {
                    "success": False,
                    "error": f"BAM文件不存在: {bam_norm}",
                }
            bam_entries.append({"sample_id": sample_id, "bam_file": bam_norm})
        if not bam_entries:
            return {"success": False, "error": "未从对齐结果中收集到任何BAM路径"}
        params_raw = (featurecounts_params or {}).copy()

        nf_params = {
            "input_bam_list": json.dumps(bam_entries, ensure_ascii=False),
            "gtf_file": gtf_file,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **params_raw,
        }

        # 资源配置：仅接受 FeatureCounts 阶段片段，规范化后注入
        from .utils_tools import normalize_resources
        normalized_fc = normalize_resources("featurecounts", {"featurecounts": resource_config or {}})
        resources_map: Dict[str, Dict[str, Any]] = {"featurecounts": normalized_fc} if normalized_fc else {}

        # M4: 参数版本化 - 使用新的接口直接传递results_dir
        versioned_params_file = None
        try:
            from .utils_tools import write_params_file
            # 直接传递results_dir，不需要创建临时state
            versioned_params_file = write_params_file(
                "featurecounts",
                nf_params,
                results_dir=str(results_dir)
            )
            logger.info(f"FeatureCounts实际运行参数版本化文件已写入: {versioned_params_file}")

        except Exception as e:
            logger.warning(f"FeatureCounts参数版本化写入失败 (不影响执行): {e}")

        # 使用版本化参数文件
        if not versioned_params_file or not versioned_params_file.exists():
            return {
                "success": False,
                "error": "版本化参数文件生成失败，无法执行FeatureCounts"
            }
        
        params_file = versioned_params_file
        logger.info(f"使用版本化参数文件运行FeatureCounts: {params_file}")

        # 定位 Nextflow 脚本
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "featurecounts.nf"
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 featurecounts.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "featurecounts.nf")],
            }

        # 执行 Nextflow
        cmd = [
            "nextflow",
            "run",
            str(nextflow_script),
            "-params-file",
            str(params_file),
            "-work-dir",
            str(work_dir),
        ]
        # 通过 -params 内联注入资源配置
        try:
            inline_params = json.dumps({"resources": resources_map}, ensure_ascii=False)
            cmd.extend(["-params", inline_params])
        except Exception as e:
            logger.warning(f"构建FeatureCounts内联资源参数失败，将不注入资源: {e}")

        logger.info("执行Nextflow FeatureCounts流水线")
        logger.info(f"参数文件: {params_file}")
        logger.info(f"工作目录: {work_dir}")
        logger.info(f"结果目录: {results_dir}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,
            cwd=tools_config.settings.project_root,
        )

        # 构建输出结构 - 聚合输出格式
        sample_count = len(bam_entries)
        fc_root = results_dir / "featurecounts"
        
        # FeatureCounts聚合输出文件
        output = {
            "counts_file": str(fc_root / "all_samples.featureCounts"),
            "summary_file": str(fc_root / "all_samples.featureCounts.summary"),
        }

        results = {
            "success": result.returncode == 0,
            "results_dir": str(results_dir),
            "sample_count": sample_count,
            "output": output
        }
        
        if results["success"]:
            logger.info(f"FeatureCounts完成：samples={sample_count} results={results_dir}")
        else:
            logger.warning(f"FeatureCounts失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        return results

    except subprocess.TimeoutExpired:
        logger.warning("FeatureCounts执行超时：30分钟")
        return {
            "success": False,
            "error": "Nextflow执行超时（30分钟）",
        }
    except Exception as e:
        logger.error(f"FeatureCounts异常：{e}")
        return {
            "success": False,
            "error": f"执行FeatureCounts流水线失败: {str(e)}",
        }


@tool
def parse_featurecounts_metrics(results_directory: str) -> Dict[str, Any]:
    """解析 FeatureCounts 定量结果为 JSON，保留样本级颗粒度

    - 读取 featurecounts/all_samples.featureCounts.summary
    - 返回每个样本的完整 Status→数值映射（不丢失任何状态）
    - 同时提供基础派生值（total_reads、assignment_rate）与总体统计

    Args:
        results_directory: FeatureCounts 结果目录（包含 featurecounts 子目录）
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {"success": False, "error": f"结果目录不存在: {results_directory}"}

        fc_dir = results_path / "featurecounts"
        if not fc_dir.exists():
            return {"success": False, "error": f"缺少特征计数目录: {fc_dir}"}

        # 汇总与计数矩阵
        summary_file = fc_dir / "all_samples.featureCounts.summary"
        counts_file = fc_dir / "all_samples.featureCounts"

        if not summary_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts汇总文件: {summary_file}"}
        if not counts_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts计数文件: {counts_file}"}




        try:
            # 使用pandas读取tab分隔的summary文件
            summary_df = pd.read_csv(summary_file, sep='	', index_col=0)
            
            if summary_df.empty:
                return {"success": False, "error": "汇总文件为空"}


            # 获取样本ID列表（列名）和状态列表（索引）
            sample_ids = summary_df.columns.tolist()
            statuses_in_file = summary_df.index.tolist()
            
            # 使用pandas直接计算统计数据
            # 添加Total行（每个样本的总reads）
            summary_df.loc['Total'] = summary_df.sum(axis=0)
            
            # 添加assignment_rate行（分配率）
            if 'Assigned' in summary_df.index:
                summary_df.loc['Assignment_Rate'] = summary_df.loc['Assigned'] / summary_df.loc['Total']
            else:
                summary_df.loc['Assignment_Rate'] = 0
            
            # 添加Overall_Total列（每个状态的总计）
            summary_df['Overall_Total'] = summary_df.drop(['Total', 'Assignment_Rate'], axis=1).sum(axis=1)
            
            # 添加Overall_Assignment_Rate行
            if 'Assigned' in summary_df.index:
                overall_total = summary_df.loc['Total', 'Overall_Total']
                overall_assigned = summary_df.loc['Assigned', 'Overall_Total'] 
                summary_df.loc['Overall_Assignment_Rate', 'Overall_Total'] = (overall_assigned / overall_total) if overall_total > 0 else 0.0


        except Exception as e:
            return {"success": False, "error": f"解析汇总文件失败: {str(e)}"}

        # 基因数（行数-表头-注释）
        gene_count = 0
        try:
            with open(counts_file, "r", encoding="utf-8") as f:
                for line in f:
                    if not line.startswith('#'):
                        gene_count += 1
                gene_count -= 1  # 减去标题行
        except Exception:
            gene_count = 0

        # 构建标准列表格式的sample_metrics
        sample_metrics = []
        for sample_id in sample_ids:
            if sample_id in summary_df.columns:
                sample_data = summary_df[sample_id].to_dict()
                sample_metrics.append({
                    "sample_id": sample_id,
                    "total_reads": sample_data.get("Total", 0),
                    "assigned_reads": sample_data.get("Assigned", 0),
                    "assignment_rate": sample_data.get("Assignment_Rate", 0.0),
                    "unassigned_nofeatures": sample_data.get("Unassigned_NoFeatures", 0),
                    "unassigned_ambiguity": sample_data.get("Unassigned_Ambiguity", 0),
                    "unassigned_mappingquality": sample_data.get("Unassigned_MappingQuality", 0),
                    "all_stats": sample_data  # 完整的统计数据
                })
        
        result = {
            "success": True,
            "results_directory": results_directory,
            "summary_file": str(summary_file),
            "counts_file": str(counts_file),
            "matrix_path": str(counts_file),
            "sample_count": len(sample_ids),
            "gene_count": gene_count,
            "statuses": statuses_in_file,
            "sample_metrics": sample_metrics,  # 标准列表格式
            "summary_data": summary_df.to_dict()  # 保留完整的pandas数据供其他用途
        }
        
        return result

    except Exception as e:
        return {
            "success": False,
            "error": f"解析FeatureCounts结果失败: {str(e)}"
        }
