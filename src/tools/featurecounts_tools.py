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
    genome_info: Dict[str, Any],
    results_timestamp: Optional[str] = None,
    base_results_dir: Optional[str] = None,
    hisat2_results: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """执行Nextflow FeatureCounts定量流程

    Args:
        featurecounts_params: FeatureCounts参数字典（可包含 -T/-s/-p/-M 等风格键）
        star_results: STAR节点结果（可为空），需包含 per_sample_outputs 中的 BAM 路径
        genome_info: 基因组信息（需提供 GTF 路径，如 gtf_path）
        results_timestamp: 可选的时间戳，优先用于结果目录
        base_results_dir: 可选的基底结果目录（来自Detect节点）
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

        # 环境检查：允许本地与容器环境，路径规范化在下方处理

        # 解析并归一化 GTF 注释文件路径（容器内）
        gtf_file_raw = (
            genome_info.get("gtf_path")
            or genome_info.get("gtf")
            or genome_info.get("annotation_gtf")
            or ""
        )
        if not gtf_file_raw:
            return {"success": False, "error": "genome_info 未提供 GTF 注释文件路径 (gtf_path)"}
        gtf_file = gtf_file_raw
        if not Path(gtf_file).exists():
            return {
                "success": False,
                "error": f"GTF文件不存在: {gtf_file}",
            }

        # 运行根目录（results_dir）：复用比对步骤的 results_dir，保持同一运行根目录
        timestamp = results_timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        run_root = Path(align_results.get("results_dir") or base_results_dir or tools_config.results_dir / f"{timestamp}")
        results_dir = run_root
        # 统一 Nextflow 工作目录到 /data/work，使用运行ID区分
        run_id = results_dir.name or timestamp
        work_dir = tools_config.settings.data_dir / "work" / f"featurecounts_{run_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)
        (results_dir / "featurecounts").mkdir(parents=True, exist_ok=True)
        logger.info(f"FeatureCounts启动：bam={len(per_sample)} results={results_dir} work={work_dir}")

        # 构建 Nextflow 参数（与 featurecounts.nf 对齐）
        # 将 STAR 输出的 BAM 列表转换为JSON字符串（Nextflow 端会 echo 后再解析）
        bam_entries = []
        for item in per_sample:
            sid = item.get("sample_id") or "sample"
            bam = item.get("aligned_bam") or item.get("bam")
            if not bam:
                continue
            bam_norm = bam
            if not Path(bam_norm).exists():
                return {
                    "success": False,
                    "error": f"BAM文件不存在: {bam_norm}",
                }
            bam_entries.append({"sample_id": sid, "bam_file": bam_norm})
        if not bam_entries:
            return {"success": False, "error": "未从STAR结果中收集到任何BAM路径"}

        # 参数映射：Python风格/短旗标 → Nextflow params 名称
        mapped: Dict[str, Any] = {}
        p = featurecounts_params or {}

        def pick_bool(key: str) -> Optional[bool]:
            v = p.get(key)
            if isinstance(v, bool):
                return v
            return None

        def pick_int(key: str) -> Optional[int]:
            v = p.get(key)
            try:
                return int(v) if v is not None else None
            except Exception:
                return None

        # 线程/链特异性/特征/属性/质量阈
        mapped["threads"] = pick_int("-T") or p.get("threads") or 4
        mapped["strand_specificity"] = pick_int("-s") or p.get("strand_specificity") or 0
        mapped["feature_type"] = p.get("-t") or p.get("feature_type") or "exon"
        mapped["attribute_type"] = p.get("-g") or p.get("attribute_type") or "gene_id"
        mapped["min_mapping_quality"] = pick_int("-Q") or p.get("min_mapping_quality") or 10

        # 布尔开关 - 修改count_reads_pairs默认值为false
        mapped["count_reads_pairs"] = pick_bool("-p") if pick_bool("-p") is not None else (p.get("count_reads_pairs") if isinstance(p.get("count_reads_pairs"), bool) else False)
        mapped["count_multi_mapping_reads"] = pick_bool("-M") if pick_bool("-M") is not None else bool(p.get("count_multi_mapping_reads", False))
        mapped["ignore_duplicates"] = bool(p.get("--ignoreDup", False) or p.get("ignore_duplicates", False))
        mapped["require_both_ends_mapped"] = bool(p.get("-B", False) or p.get("require_both_ends_mapped", False))
        mapped["exclude_chimeric"] = bool(p.get("-C", False) or p.get("exclude_chimeric", False))

        # 组装 Nextflow 参数文件
        nf_params = {
            "input_bam_list": json.dumps(bam_entries, ensure_ascii=False),
            "gtf_file": gtf_file,
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            **mapped,
        }

        # M4: 参数版本化 - 将实际用于运行的nf_params写入版本化文件
        versioned_params_file = None
        try:
            # 构建临时 state 用于参数版本化
            from ..state import AgentState
            from .utils_tools import write_params_file
            temp_state = AgentState()
            temp_state.results_dir = str(results_dir)
            temp_state.results_timestamp = results_timestamp or ""
            
            # 写入版本化文件 - 使用实际运行的nf_params
            versioned_params_file = write_params_file("featurecounts", nf_params, temp_state)
            logger.info(f"FeatureCounts实际运行参数版本化文件已写入: {versioned_params_file}")
            
        except Exception as e:
            logger.warning(f"FeatureCounts参数版本化写入失败 (不影响执行): {e}")

        # 使用版本化文件作为参数文件，如果失败则降级到临时文件
        if versioned_params_file and versioned_params_file.exists():
            params_file = versioned_params_file
            logger.info(f"使用版本化参数文件运行FeatureCounts: {params_file}")
        else:
            # 降级方案：生成临时参数文件
            logger.warning("版本化参数文件不可用，使用临时参数文件")
            fc_dir = results_dir / "featurecounts"
            fc_dir.mkdir(parents=True, exist_ok=True)
            params_file = fc_dir / "featurecounts_params.json"
            with open(params_file, "w", encoding="utf-8") as f:
                json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 定位 Nextflow 脚本
        nextflow_script = Path('/src/nextflow/featurecounts.nf')
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 featurecounts.nf",
                "searched": ["/src/nextflow/featurecounts.nf"],
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

        # 构建输出结构 - 适配新的批量输出格式
        sample_count = len(bam_entries)
        per_sample_outputs: List[Dict[str, Any]] = []
        fc_root = results_dir / "featurecounts"
        
        # 新的featurecounts.nf脚本生成批量文件，不再有每个样本的单独目录
        # 主要输出文件：
        # - all_samples.featureCounts (完整计数矩阵)
        # - all_samples.featureCounts.summary (统计摘要)
        # - merged_counts_matrix.txt (兼容格式的矩阵)
        
        # 为兼容性生成per_sample_outputs结构，指向批量文件
        for entry in bam_entries:
            sid = entry["sample_id"]
            sample_output = {
                "sample_id": sid,
                "counts_file": str(fc_root / "all_samples.featureCounts"),  # 指向批量文件
                "summary_file": str(fc_root / "all_samples.featureCounts.summary"),  # 指向批量文件
            }
            per_sample_outputs.append(sample_output)

        payload = {
            "success": result.returncode == 0,
            "message": f"FeatureCounts定量完成，处理了{sample_count}个样本" if result.returncode == 0 else "FeatureCounts执行失败",
            "results_dir": str(results_dir),
            "work_dir": str(work_dir),
            "params_file": str(params_file),
            "sample_count": sample_count,
            "per_sample_outputs": per_sample_outputs,
            "matrix_path": str(fc_root / "merged_counts_matrix.txt"),
            "nextflow_params": nf_params,
        }
        try:
            if get_tools_config().settings.debug_mode:
                payload.update({
                    "stdout": result.stdout,
                    "stderr": result.stderr,
                    "cmd": " ".join(cmd),
                })
            else:
                # 非调试模式去掉 nextflow_params 以减小负载
                payload.pop("nextflow_params", None)
        except Exception:
            pass
        if payload["success"]:
            logger.info(f"FeatureCounts完成：samples={sample_count} results={results_dir}")
        else:
            logger.warning(f"FeatureCounts失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        return payload

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
    """解析FeatureCounts定量结果，输出样本级与总体指标
    
    Args:
        results_directory: FeatureCounts结果目录（包含 featurecounts 子目录）
    
    Returns:
        解析后的指标（assignment rates、未分配原因等）
    """
    try:
        results_path = Path(results_directory)
        if not results_path.exists():
            return {"success": False, "error": f"结果目录不存在: {results_directory}"}
        
        fc_dir = results_path / "featurecounts"
        if not fc_dir.exists():
            return {"success": False, "error": f"缺少特征计数目录: {fc_dir}"}
        
        # 查找批量输出的汇总文件
        summary_file = fc_dir / "all_samples.featureCounts.summary"
        counts_file = fc_dir / "all_samples.featureCounts"
        
        if not summary_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts汇总文件: {summary_file}"}
        
        if not counts_file.exists():
            return {"success": False, "error": f"未找到FeatureCounts计数文件: {counts_file}"}
        
        # 解析批量汇总文件
        sample_metrics: List[Dict[str, Any]] = []
        totals = {
            "assigned": 0,
            "nofeatures": 0,
            "multimapping": 0,
            "ambiguous": 0,
            "mappingquality": 0,
            "other": 0,
            "total": 0,
        }
        
        try:
            # 规范化样本ID的内部工具：
            # - 兼容列名为BAM文件路径/文件名/带STAR后缀的多种情况
            # - 目标：与样本ID（如 SRRxxxx、样本目录名）对齐
            def _normalize_sample_id(name: str) -> str:
                s = str(name or "").strip()
                if not s:
                    return s
                # 去除可能的路径前缀（同时支持 / 与 \\ 分隔符）
                if "/" in s:
                    s = s.split("/")[-1]
                if "\\" in s:
                    s = s.split("\\")[-1]
                # 去除常见扩展名
                for ext in [".bam", ".cram", ".sam", ".txt"]:
                    if s.endswith(ext):
                        s = s[: -len(ext)]
                        break
                # 去除常见后缀（STAR/HISAT2 的命名后缀，点/下划线变体）
                star_suffixes = [
                    ".Aligned.sortedByCoord.out",
                    ".Aligned.out",
                    ".Aligned",
                    "_Aligned.sortedByCoord.out",
                    "_Aligned.out",
                    "_Aligned",
                    ".hisat2",  # 来自 HISAT2 的常见后缀（在移除 .bam 后可能残留）
                    "_hisat2",
                ]
                for suf in star_suffixes:
                    if s.endswith(suf):
                        s = s[: -len(suf)]
                        break
                return s

            with open(summary_file, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
                
                if len(lines) < 2:
                    return {"success": False, "error": "汇总文件格式错误"}
                
                # 解析标题行获取样本名称（用于回退）
                header = lines[0].strip().split('\t')
                if len(header) < 2:
                    return {"success": False, "error": "汇总文件标题行格式错误"}

                # 优先从参数文件读取样本ID顺序（与执行输入一致）
                sample_ids: List[str] = []
                try:
                    params_path = results_path / "featurecounts" / "featurecounts_params.json"
                    if params_path.exists():
                        with open(params_path, "r", encoding="utf-8") as pf:
                            pf_json = json.load(pf)
                        input_bam_list = pf_json.get("input_bam_list")
                        if isinstance(input_bam_list, str):
                            input_bam_list = json.loads(input_bam_list)
                        if isinstance(input_bam_list, list):
                            for ent in input_bam_list:
                                sid = ent.get("sample_id") or _normalize_sample_id(ent.get("bam_file", ""))
                                sample_ids.append(_normalize_sample_id(sid))
                except Exception:
                    # 若读取或解析失败，忽略并回退到header
                    sample_ids = []

                # 校验样本数是否与汇总列数一致；否则回退到header列名
                if not sample_ids or len(sample_ids) != (len(header) - 1):
                    sample_names = header[1:]  # 第一列是Status，后面是样本名（通常为输入BAM的文件名）
                    sample_ids = [_normalize_sample_id(nm) for nm in sample_names]

                # 初始化每个样本的指标
                for sid in sample_ids:
                    
                    sample_metrics.append({
                        "sample_id": sid,
                        "assigned": 0,
                        "unassigned_unmapped": 0,
                        "unassigned_mappingquality": 0,
                        "unassigned_nofeatures": 0,
                        "unassigned_ambiguity": 0,
                        "total_reads": 0
                    })
                
                # 解析每一行统计数据
                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        continue
                        
                    status = parts[0]
                    values = [int(v) for v in parts[1:]]
                    
                    # 更新每个样本的指标
                    for i, value in enumerate(values):
                        if i < len(sample_metrics):
                            if status == "Assigned":
                                sample_metrics[i]["assigned"] = value
                                totals["assigned"] += value
                            elif status == "Unassigned_Unmapped":
                                sample_metrics[i]["unassigned_unmapped"] = value
                            elif status == "Unassigned_MappingQuality":
                                sample_metrics[i]["unassigned_mappingquality"] = value
                                totals["mappingquality"] += value
                            elif status == "Unassigned_NoFeatures":
                                sample_metrics[i]["unassigned_nofeatures"] = value
                                totals["nofeatures"] += value
                            elif status == "Unassigned_Ambiguity":
                                sample_metrics[i]["unassigned_ambiguity"] = value
                                totals["ambiguous"] += value
                
                # 计算每个样本的总读数和分配率
                for sample_metric in sample_metrics:
                    sample_metric["total_reads"] = (
                        sample_metric["assigned"] +
                        sample_metric["unassigned_unmapped"] +
                        sample_metric["unassigned_mappingquality"] +
                        sample_metric["unassigned_nofeatures"] +
                        sample_metric["unassigned_ambiguity"]
                    )
                    
                    if sample_metric["total_reads"] > 0:
                        sample_metric["assignment_rate"] = round(
                            sample_metric["assigned"] / sample_metric["total_reads"], 4
                        )
                    else:
                        sample_metric["assignment_rate"] = 0.0
                    
                    totals["total"] += sample_metric["total_reads"]
        
        except Exception as e:
            return {"success": False, "error": f"解析汇总文件失败: {str(e)}"}
        
        # 计算总体统计
        total_assignment_rate = totals["assigned"] / totals["total"] if totals["total"] > 0 else 0.0
        
        # 读取计数矩阵获取基因数量
        gene_count = 0
        try:
            with open(counts_file, "r", encoding="utf-8") as f:
                # 跳过注释行
                for line in f:
                    if not line.startswith('#'):
                        gene_count += 1
                gene_count -= 1  # 减去标题行
        except Exception:
            gene_count = 0
        
        # 重要文件路径（若存在合并矩阵，则优先提供）
        matrix_file = fc_dir / "merged_counts_matrix.txt"

        return {
            "success": True,
            "results_directory": results_directory,
            "summary_file": str(summary_file),
            "counts_file": str(counts_file),
            "matrix_path": str(matrix_file) if matrix_file.exists() else str(counts_file),
            "sample_count": len(sample_metrics),
            "gene_count": gene_count,
            "sample_metrics": sample_metrics,
            "overall_statistics": {
                "total_reads": totals["total"],
                "total_assigned": totals["assigned"],
                "total_unassigned_nofeatures": totals["nofeatures"],
                "total_unassigned_ambiguity": totals["ambiguous"],
                "total_unassigned_mappingquality": totals["mappingquality"],
                "overall_assignment_rate": round(total_assignment_rate, 4)
            }
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": f"解析FeatureCounts结果失败: {str(e)}"
        }