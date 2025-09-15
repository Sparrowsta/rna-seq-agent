"""
Analysis节点  - 基于设计文档的综合分析实现

功能：
- 解析 FastP/STAR/FeatureCounts 三步结果
- 计算样本健康度和总体结论
- 生成结构化 JSON 报告和 Markdown 摘要
- 集成 LLM 智能分析
- 支持历史归档和报告落盘
"""

import json
import asyncio
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from ..state import AgentState, LLMAnalysisModel
from ..prompts import ANALYSIS_LLM_SYSTEM_PROMPT
from ..tools import (
    parse_fastp_results, 
    parse_star_metrics, 
    parse_featurecounts_metrics,
    write_analysis_markdown
)
from ..core import get_shared_llm
from ..logging_bootstrap import get_logger, log_llm_preview

logger = get_logger("rna.analysis")


async def analysis_node(state: AgentState) -> Dict[str, Any]:
    """
    Analysis节点  - 综合分析实现
    
    执行流程：
    1. 前置校验 - 检查所有步骤结果
    2. 解析三步指标 - 调用解析器工具
    3. 样本ID归一化和对齐
    4. 指标合并与健康度评估
    5. LLM智能总结（可选）
    6. 报告落盘 - JSON和Markdown
    7. 状态回填与清理
    """
    logger.info("综合分析节点开始执行")
    
    try:
        # 1. 前置校验
        validation_result = _validate_input_results(state)
        if not validation_result["success"]:
            return _create_error_response(validation_result["error"])
        
        results_dir = validation_result["results_dir"]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        logger.info(f"分析结果目录: {results_dir}")
        
        # 2. 解析三步指标
        logger.info("解析 FastP/STAR/FeatureCounts 结果...")
        parsing_result = _parse_pipeline_metrics(results_dir)
        
        if not parsing_result["success"]:
            return _create_error_response(f"指标解析失败: {parsing_result['error']}")
        
        fastp_data = parsing_result["fastp"]
        star_data = parsing_result["star"] 
        featurecounts_data = parsing_result["featurecounts"]
        parsing_status = parsing_result.get("parsing_status", {})
        parsing_errors = parsing_result.get("parsing_errors", {})
        
        # 3. 样本ID归一化和对齐
        logger.debug("样本ID归一化和指标对齐...")
        alignment_result = _align_sample_metrics(
            state.nextflow_config.get("sample_groups", []),
            fastp_data, star_data, featurecounts_data
        )
        
        # 4. 指标合并与健康度评估
        logger.debug("计算样本健康度和总体结论...")
        assessment_result = _assess_sample_health(alignment_result)
        
        # 5. 构建基础报告结构
        base_report = _build_base_report(
            state, results_dir, timestamp, 
            fastp_data, star_data, featurecounts_data,
            alignment_result, assessment_result
        )
        
        # 添加解析状态到报告中，增强可观测性
        if parsing_errors:
            base_report.setdefault("files", {})["parsing_errors"] = parsing_errors
            failed_steps = [step for step, success in parsing_status.items() if not success]
            if failed_steps:
                base_report.setdefault("summary", {}).setdefault("key_findings", []).append(
                    f"⚠️ 解析失败的步骤: {', '.join(failed_steps)}")
                logger.warning(f"部分解析失败: {failed_steps}")
        
        # 6. LLM智能总结（可选但推荐）
        logger.info("执行LLM智能分析...")
        llm_result = await _execute_llm_analysis(base_report)
        
        # 将LLM结果合并到报告中，并增强可观测性
        if llm_result["success"]:
            base_report["llm"] = llm_result["analysis"]
            logger.info("LLM智能分析完成")
        else:
            base_report["llm_error"] = llm_result["error"]
            # 在key_findings中添加LLM降级提示，提升用户可见性
            base_report.setdefault("summary", {}).setdefault("key_findings", []).insert(0, 
                f"🤖 LLM智能分析不可用（降级）：{llm_result['error']}")
            logger.warning(f"LLM分析失败，使用规则分级结果: {llm_result['error']}")
        
        # 7. 报告落盘
        logger.info("生成并保存分析报告...")
        report_result = _save_reports(base_report, results_dir, timestamp)
        
        if not report_result["success"]:
            logger.error(f"报告保存失败: {report_result['error']}")
        
        # 8. 状态回填与清理
        logger.debug("更新状态并清理...")
        return _create_success_response(base_report, report_result)
        
    except Exception as e:
        logger.error(f"Analysis Node 执行异常: {str(e)}", exc_info=True)
        return _create_error_response(f"分析节点执行异常: {str(e)}")


async def _invoke_llm_langgraph(structured_llm, messages: List[Dict[str, Any]]):
    """以 LangGraph 风格优先的方式异步调用 LLM（结构化输出）。

    优先尝试传入 {"messages": [...]}；
    若模型不接受该输入格式，则回退为纯列表 [...]
    """
    try:
        return await structured_llm.ainvoke({"messages": messages})
    except Exception:
        # 回退为直接传列表
        return await structured_llm.ainvoke(messages)


def _validate_input_results(state: AgentState) -> Dict[str, Any]:
    """前置校验：检查所有步骤的结果状态"""
    fastp_results = state.fastp_results or {}
    star_results = state.star_results or {}
    featurecounts_results = state.featurecounts_results or {}
    
    # 检查至少有一个步骤成功
    any_success = (
        fastp_results.get("success") or
        star_results.get("success") or  
        featurecounts_results.get("success")
    )
    
    if not any_success:
        return {
            "success": False,
            "error": "所有流水线步骤均未成功执行，无法进行综合分析"
        }
    
    # 确定结果目录优先级：featurecounts > star > fastp
    results_dir = None
    for results in [featurecounts_results, star_results, fastp_results]:
        if results.get("success") and results.get("results_dir"):
            results_dir = results["results_dir"]
            break
    
    if not results_dir:
        return {
            "success": False,
            "error": "无法确定分析结果目录路径"
        }
    
    if not Path(results_dir).exists():
        return {
            "success": False,
            "error": f"结果目录不存在: {results_dir}"
        }
    
    return {
        "success": True,
        "results_dir": results_dir
    }


def _parse_pipeline_metrics(results_dir: str) -> Dict[str, Any]:
    """解析三个流水线步骤的指标，增强错误可观测性"""
    parsing_status = {"fastp": False, "star": False, "featurecounts": False}
    parsing_errors = {}
    
    try:
        # 调用解析器工具，分别捕获每个解析器的成功/失败状态
        logger.debug("解析FastP结果...")
        try:
            fastp_result = parse_fastp_results.invoke({"results_directory": results_dir})
            parsing_status["fastp"] = fastp_result.get("success", False)
            if not parsing_status["fastp"]:
                parsing_errors["fastp"] = fastp_result.get("error", "未知FastP解析错误")
        except Exception as e:
            fastp_result = {"success": False, "error": f"FastP解析异常: {str(e)}"}
            parsing_errors["fastp"] = str(e)
            
        logger.debug("解析STAR结果...")
        try:
            star_result = parse_star_metrics.invoke({"results_directory": results_dir})
            parsing_status["star"] = star_result.get("success", False)
            if not parsing_status["star"]:
                parsing_errors["star"] = star_result.get("error", "未知STAR解析错误")
        except Exception as e:
            star_result = {"success": False, "error": f"STAR解析异常: {str(e)}"}
            parsing_errors["star"] = str(e)
            
        logger.debug("解析FeatureCounts结果...")
        try:
            fc_result = parse_featurecounts_metrics.invoke({"results_directory": results_dir})
            parsing_status["featurecounts"] = fc_result.get("success", False)
            if not parsing_status["featurecounts"]:
                parsing_errors["featurecounts"] = fc_result.get("error", "未知FeatureCounts解析错误")
        except Exception as e:
            fc_result = {"success": False, "error": f"FeatureCounts解析异常: {str(e)}"}
            parsing_errors["featurecounts"] = str(e)
        
        # 统计解析成功情况
        success_count = sum(parsing_status.values())
        total_count = len(parsing_status)
        
        logger.info(f"解析完成: {success_count}/{total_count} 个步骤成功")
        if parsing_errors:
            logger.warning(f"解析错误: {list(parsing_errors.keys())}")
        
        return {
            "success": True,
            "fastp": fastp_result,
            "star": star_result, 
            "featurecounts": fc_result,
            "parsing_status": parsing_status,
            "parsing_errors": parsing_errors
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"调用解析器工具失败: {str(e)}",
            "parsing_status": parsing_status,
            "parsing_errors": parsing_errors
        }


def _align_sample_metrics(sample_groups: List[Dict], fastp_data: Dict, star_data: Dict, fc_data: Dict) -> Dict[str, Any]:
    """样本指标对齐（简化版）

    依赖各步骤返回的 sample_id 已一致（FeatureCounts 由参数文件提供权威样本顺序），
    直接按 sample_groups 中的 sample_id 精确合并，不再做额外规范化回退。
    """

    # 提取预期的样本ID列表（保持顺序）
    expected_sample_ids: List[str] = []
    for group in sample_groups:
        sid = group.get("sample_id") or group.get("id")
        if sid:
            expected_sample_ids.append(sid)

    # 各步骤样本指标映射
    fastp_samples = {s.get("sample_id"): s for s in fastp_data.get("sample_metrics", []) if s.get("sample_id")}
    star_samples = {s.get("sample_id"): s for s in star_data.get("sample_metrics", []) if s.get("sample_id")}
    fc_samples = {s.get("sample_id"): s for s in fc_data.get("sample_metrics", []) if s.get("sample_id")}

    # 合并
    aligned_samples: List[Dict[str, Any]] = []
    for sid in expected_sample_ids:
        sample_data = {
            "sample_id": sid,
            "fastp": fastp_samples.get(sid, {"error": "未找到FastP数据"}),
            "star": star_samples.get(sid, {"error": "未找到STAR数据"}),
            "featurecounts": fc_samples.get(sid, {"error": "未找到FeatureCounts数据"}),
            "notes": []
        }

        if "error" in sample_data["fastp"]:
            sample_data["notes"].append("FastP数据缺失")
        if "error" in sample_data["star"]:
            sample_data["notes"].append("STAR数据缺失")
        if "error" in sample_data["featurecounts"]:
            sample_data["notes"].append("FeatureCounts数据缺失")

        aligned_samples.append(sample_data)

    return {
        "samples": aligned_samples,
        "expected_count": len(expected_sample_ids),
        "fastp_available": len(fastp_samples),
        "star_available": len(star_samples),
        "featurecounts_available": len(fc_samples)
    }


def _assess_sample_health(alignment_data: Dict[str, Any]) -> Dict[str, Any]:
    """计算样本健康度和总体结论"""
    
    # 阈值定义（基于设计文档第3节）
    thresholds = {
        "fastp": {
            "q30_good": 0.85, "q30_warn": 0.7,
            "pass_rate_good": 0.8, "pass_rate_warn": 0.6
        },
        "star": {
            "mapping_good": 0.85, "mapping_warn": 0.7,
            "unique_good": 0.8, "unique_warn": 0.6,
            "mismatch_good": 0.05, "mismatch_warn": 0.08
        },
        "featurecounts": {
            "assignment_good": 0.6, "assignment_warn": 0.4
        }
    }
    
    def _evaluate_metric(value, good_thresh, warn_thresh, reverse=False):
        """评估单个指标（reverse=True表示越小越好）"""
        if value is None:
            return "UNKNOWN"
        try:
            val = float(value)
            if reverse:
                return "PASS" if val <= good_thresh else ("WARN" if val <= warn_thresh else "FAIL")
            else:
                return "PASS" if val >= good_thresh else ("WARN" if val >= warn_thresh else "FAIL")
        except (ValueError, TypeError):
            return "UNKNOWN"
    
    assessed_samples = []
    health_counts = {"PASS": 0, "WARN": 0, "FAIL": 0, "UNKNOWN": 0}
    
    for sample in alignment_data["samples"]:
        sample_id = sample["sample_id"]
        fastp_metrics = sample.get("fastp", {})
        star_metrics = sample.get("star", {})  
        fc_metrics = sample.get("featurecounts", {})
        
        # 评估各步骤
        evaluations = []
        
        # FastP评估
        if "error" not in fastp_metrics:
            q30_status = _evaluate_metric(
                fastp_metrics.get("q30_after"), 
                thresholds["fastp"]["q30_good"], 
                thresholds["fastp"]["q30_warn"]
            )
            pass_rate_status = _evaluate_metric(
                fastp_metrics.get("read_pass_rate"),
                thresholds["fastp"]["pass_rate_good"],
                thresholds["fastp"]["pass_rate_warn"] 
            )
            fastp_status = "FAIL" if "FAIL" in [q30_status, pass_rate_status] else (
                "WARN" if "WARN" in [q30_status, pass_rate_status] else "PASS"
            )
            evaluations.append(fastp_status)
        
        # STAR评估
        if "error" not in star_metrics:
            mapping_status = _evaluate_metric(
                star_metrics.get("mapping_rate"),
                thresholds["star"]["mapping_good"],
                thresholds["star"]["mapping_warn"]
            )
            unique_status = _evaluate_metric(
                star_metrics.get("unique_mapping_rate"),
                thresholds["star"]["unique_good"], 
                thresholds["star"]["unique_warn"]
            )
            mismatch_status = _evaluate_metric(
                star_metrics.get("mismatch_rate"),
                thresholds["star"]["mismatch_good"],
                thresholds["star"]["mismatch_warn"],
                reverse=True
            )
            star_status = "FAIL" if "FAIL" in [mapping_status, unique_status, mismatch_status] else (
                "WARN" if "WARN" in [mapping_status, unique_status, mismatch_status] else "PASS"
            )
            evaluations.append(star_status)
        
        # FeatureCounts评估
        if "error" not in fc_metrics:
            assignment_status = _evaluate_metric(
                fc_metrics.get("assignment_rate"),
                thresholds["featurecounts"]["assignment_good"],
                thresholds["featurecounts"]["assignment_warn"]
            )
            evaluations.append(assignment_status)
        
        # 综合健康度：取最差级别
        if not evaluations:
            overall_health = "UNKNOWN"
        else:
            if "FAIL" in evaluations:
                overall_health = "FAIL" 
            elif "WARN" in evaluations:
                overall_health = "WARN"
            else:
                overall_health = "PASS"
        
        health_counts[overall_health] += 1
        
        assessed_sample = {
            "sample_id": sample_id,
            "health": overall_health,
            "fastp": fastp_metrics,
            "star": star_metrics,
            "featurecounts": fc_metrics,
            "notes": sample.get("notes", [])
        }
        assessed_samples.append(assessed_sample)
    
    # 总体结论
    total_samples = len(assessed_samples)
    if health_counts["FAIL"] > 0:
        overall_status = "FAIL" if health_counts["FAIL"] / total_samples > 0.5 else "WARN"
    elif health_counts["WARN"] > 0:
        overall_status = "WARN"
    else:
        overall_status = "PASS"
    
    return {
        "samples": assessed_samples,
        "summary": {
            "status": overall_status,
            "samples": {
                "total": total_samples,
                "pass": health_counts["PASS"],
                "warn": health_counts["WARN"], 
                "fail": health_counts["FAIL"],
                "unknown": health_counts["UNKNOWN"]
            }
        }
    }


def _build_base_report(state: AgentState, results_dir: str, timestamp: str, 
                      fastp_data: Dict, star_data: Dict, fc_data: Dict,
                      alignment_result: Dict, assessment_result: Dict) -> Dict[str, Any]:
    """构建基础报告结构"""
    
    nextflow_config = state.nextflow_config or {}
    
    return {
        "pipeline": {
            "steps": ["fastp", "star", "featurecounts"],
            "species": nextflow_config.get("species", "unknown"),
            "genome_version": nextflow_config.get("genome_version", "unknown")
        },
        "context": {
            "results_dir": results_dir,
            "timestamp": timestamp,
            "sample_count": len(alignment_result.get("samples", []))
        },
        "metrics": {
            "fastp": {
                "overall": fastp_data.get("overall_statistics", {}),
                "samples": fastp_data.get("sample_metrics", [])
            },
            "star": {
                "overall": star_data.get("overall_statistics", {}),
                "samples": star_data.get("sample_metrics", [])
            },
            "featurecounts": {
                "overall": fc_data.get("overall_statistics", {}),
                "samples": fc_data.get("sample_metrics", [])
            }
        },
        "per_sample": assessment_result["samples"],
        "summary": assessment_result["summary"],
        "files": {
            "matrix_path": fc_data.get("matrix_path", "")
        },
        "recommendations": _generate_basic_recommendations(assessment_result)
    }


def _generate_basic_recommendations(assessment_result: Dict) -> List[Dict[str, str]]:
    """生成基础建议（规则基础，LLM可补充）"""
    recommendations = []
    
    summary = assessment_result.get("summary", {})
    samples_info = summary.get("samples", {})
    
    # 基于健康度统计的建议
    if samples_info.get("fail", 0) > 0:
        recommendations.append({
            "type": "action",
            "title": "检查失败样本",
            "detail": f"有{samples_info['fail']}个样本质量评估为FAIL，建议检查原始数据质量、流水线参数设置"
        })
    
    if samples_info.get("warn", 0) > 0:
        recommendations.append({
            "type": "action", 
            "title": "关注警告样本",
            "detail": f"有{samples_info['warn']}个样本存在质量警告，建议优化参数或排除异常样本"
        })
    
    # 通用后续分析建议
    if samples_info.get("pass", 0) > 0:
        recommendations.append({
            "type": "next",
            "title": "差异表达分析", 
            "detail": "可使用DESeq2或edgeR等工具进行差异表达分析"
        })
        
        recommendations.append({
            "type": "next",
            "title": "功能富集分析",
            "detail": "建议进行GO enrichment和KEGG pathway分析"
        })
    
    return recommendations


async def _execute_llm_analysis(base_report: Dict[str, Any]) -> Dict[str, Any]:
    """执行LLM智能分析，包含退避和降级策略"""
    
    # 预定义变量避免UnboundLocalError
    llm = None
    structured_llm = None 
    system_prompt = ""
    user_message = ""
    
    try:
        llm = get_shared_llm()
        
        # 构建系统提示
        system_prompt = ANALYSIS_LLM_SYSTEM_PROMPT

        # 构建用户消息（异常样本优先采样）
        pipeline_info = base_report.get("pipeline", {})
        context = base_report.get("context", {})
        summary = base_report.get("summary", {})
        per_sample = base_report.get("per_sample", [])
        
        # 异常样本优先采样策略
        sorted_samples = sorted(per_sample, key=lambda x: {
            "FAIL": 0, "WARN": 1, "PASS": 2, "UNKNOWN": 3
        }.get(x.get("health", "UNKNOWN"), 3))
        sampled_per_sample = sorted_samples[:10]  # 限制样本数量避免过长
        
        user_message = f"""
分析基本信息：
- 流程步骤：{' → '.join(pipeline_info.get('steps', []))}
- 物种：{pipeline_info.get('species')}，基因组版本：{pipeline_info.get('genome_version')}
- 样本数量：{context.get('sample_count')}
- 分析目录：{context.get('results_dir')}

总体结论：
- 状态：{summary.get('status')}
- 样本统计：{summary.get('samples')}

前{len(sampled_per_sample)}个样本详情（异常样本优先）：
{json.dumps(sampled_per_sample, ensure_ascii=False, indent=2)}

请提供：
1. 面向用户的详细的总体评估
2. 关键发现列表（包含具体数据）
3. 每个有问题样本的具体问题和严重性
4. 分类建议（操作建议、后续分析方向） 
5. 潜在风险提示
"""

        # 调用LLM（带超时）
        structured_llm = llm.with_structured_output(LLMAnalysisModel)
        
        logger.info("调用LLM进行智能分析...")
        msgs = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_message}
        ]
        llm_response = await _invoke_llm_langgraph(structured_llm, msgs)
        try:
            log_llm_preview(logger, "analysis", llm_response)
        except Exception:
            pass
        
        return {
            "success": True,
            "analysis": dict(llm_response)
        }
        
    except Exception as e:
        # 实现退避策略
        error_msg = str(e)
        
        # 更精确的错误码识别
        if "429" in error_msg:  # 速率限制
            logger.warning("遇到速率限制，等待20秒后重试...")
            await asyncio.sleep(20)
            try:
                # 重新构造LLM和消息，避免引用未定义变量
                if llm is None or structured_llm is None:
                    llm = get_shared_llm()
                    structured_llm = llm.with_structured_output(LLMAnalysisModel)
                if not system_prompt:
                    system_prompt = ANALYSIS_LLM_SYSTEM_PROMPT
                if not user_message:
                    # 重新构造user_message（简化版避免重复逻辑）
                    user_message = f"简化分析请求 - 状态：{base_report.get('summary', {}).get('status', 'UNKNOWN')}"
                
                msgs = [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_message}
                ]
                llm_response = await _invoke_llm_langgraph(structured_llm, msgs)
                try:
                    log_llm_preview(logger, "analysis.retry429", llm_response)
                except Exception:
                    pass
                return {
                    "success": True,
                    "analysis": dict(llm_response)
                }
            except Exception as retry_e:
                return {
                    "success": False,
                    "error": f"LLM分析失败（速率限制重试后）：{str(retry_e)}"
                }
        
        # 更精确的5xx服务器错误识别
        elif any(code in error_msg for code in ["500", "502", "503", "504", "timeout"]):  # 服务器错误或超时
            logger.warning("遇到服务器问题，等待2秒后重试...")
            await asyncio.sleep(2)
            try:
                # 重新构造LLM和消息
                if llm is None or structured_llm is None:
                    llm = get_shared_llm()
                    structured_llm = llm.with_structured_output(LLMAnalysisModel)
                if not system_prompt:
                    system_prompt = ANALYSIS_LLM_SYSTEM_PROMPT
                if not user_message:
                    user_message = f"简化分析请求 - 状态：{base_report.get('summary', {}).get('status', 'UNKNOWN')}"
                
                msgs = [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_message}
                ]
                llm_response = await _invoke_llm_langgraph(structured_llm, msgs)
                try:
                    log_llm_preview(logger, "analysis.retry5xx", llm_response)
                except Exception:
                    pass
                return {
                    "success": True,
                    "analysis": dict(llm_response)
                }
            except Exception as retry_e:
                return {
                    "success": False,
                    "error": f"LLM分析失败（服务器错误重试后）：{str(retry_e)}"
                }
        
        return {
            "success": False,
            "error": f"LLM分析失败：{error_msg}"
        }


def _save_reports(report_data: Dict[str, Any], results_dir: str, timestamp: str) -> Dict[str, Any]:
    """保存JSON报告和Markdown摘要"""
    
    try:
        # 创建报告目录
        report_dir = Path(results_dir) / "reports" / timestamp
        report_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存JSON报告
        json_path = report_dir / "analysis_report.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=2, ensure_ascii=False)
        
        # 更新报告中的文件路径
        report_data.setdefault("files", {})["report_json"] = str(json_path)
        
        # 调用Markdown生成工具
        md_result = write_analysis_markdown.invoke({
            "analysis_report": report_data,
            "output_dir": str(report_dir),
            "filename": "analysis_summary.md",
            "append_llm_section": True
        })
        
        if md_result.get("success"):
            report_data["files"]["report_md"] = md_result["path"]
        
        return {
            "success": True,
            "json_path": str(json_path),
            "md_path": md_result.get("path") if md_result.get("success") else None,
            "report_dir": str(report_dir)
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": f"报告保存失败：{str(e)}"
        }


def _create_success_response(report_data: Dict[str, Any], save_result: Dict[str, Any]) -> Dict[str, Any]:
    """创建成功响应"""
    
    # 生成用户友好的摘要
    summary = report_data.get("summary", {})
    
    status = summary.get("status", "UNKNOWN")
    status_emoji = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}.get(status, "❓")
    
    samples_info = summary.get("samples", {})
    
    user_response = f"""
🎉 RNA-seq综合分析完成！

{status_emoji} 总体结论: {status}

📊 样本统计:
- 总计：{samples_info.get('total', 0)} 个样本
- 通过：{samples_info.get('pass', 0)} 个 ✅
- 警告：{samples_info.get('warn', 0)} 个 ⚠️  
- 失败：{samples_info.get('fail', 0)} 个 ❌

📁 分析报告:
- JSON详细报告: {save_result.get('json_path', '生成失败')}
- Markdown摘要: {save_result.get('md_path', '生成失败')}

💡 后续建议: 查看详细报告了解具体指标和建议
"""

    # 提取关键统计信息
    workflow_statistics = {
        "total_samples": samples_info.get("total", 0),
        "pass_samples": samples_info.get("pass", 0), 
        "warn_samples": samples_info.get("warn", 0),
        "fail_samples": samples_info.get("fail", 0),
        "overall_status": status
    }
    
    return {
        "success": True,
        "status": "analysis_completed",
        "response": user_response.strip(),
        "analysis_report": report_data,
        "analysis_report_path": save_result.get("json_path"),
        "workflow_statistics": workflow_statistics,
        
        # 清空执行进度状态
        "current_step": "",
        "completed_steps": [],
        "execution_mode": "single",
        
        # 清空各节点的结果状态（但保留引用）
        "fastp_results": {},
        "star_results": {},
        "featurecounts_results": {},
        
        # 清空优化相关状态
        "fastp_optimization_suggestions": "",
        "star_optimization_suggestions": "",
        "featurecounts_optimization_suggestions": ""
    }


def _create_error_response(error_message: str) -> Dict[str, Any]:
    """创建错误响应"""
    return {
        "success": False,
        "status": "analysis_failed",
        "response": f"❌ 综合分析失败\n\n错误详情: {error_message}\n\n建议：请检查上游流水线执行状态，确保至少有一个步骤成功完成。",
        "analysis_results": {
            "success": False,
            "status": "failed", 
            "error": error_message
        }
    }
