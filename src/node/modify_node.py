"""
Modify Node - 智能配置修改节点
负责解析用户修改需求并更新所有相关配置参数
"""
from typing import Dict, Any, List
from pathlib import Path
from datetime import datetime
from pydantic import BaseModel, Field
from ..state import AgentState
from ..logging_bootstrap import get_logger, log_llm_preview
from ..core import get_shared_llm
from ..prompts import MODIFY_NODE_PROMPT
import json


class ModifyRequest(BaseModel):
    """修改请求的结构化输出"""
    nextflow_changes: Dict[str, Any] = Field(
        default={}, 
        description="Nextflow配置修改：species, genome_version, qc_tool, align_tool, quant_tool等"
    )
    resource_changes: Dict[str, Dict[str, Any]] = Field(
        default={}, 
        description="资源配置修改：各进程的cpus和memory设置"
    )
    fastp_changes: Dict[str, Any] = Field(
        default={}, 
        description="FastP参数修改：quality_threshold, length_required, adapter_trimming等"
    )
    star_changes: Dict[str, Any] = Field(
        default={},
        description="STAR参数修改：outFilterMultimapNmax, twopassMode, quantMode等STAR特有参数。当用户明确提到STAR或这些参数时必须使用此字段！"
    )
    hisat2_changes: Dict[str, Any] = Field(
        default={},
        description="HISAT2参数修改：--mp, --rdg, --rfg, --score-min等HISAT2特有参数。当用户明确提到HISAT2或这些参数时必须使用此字段！"
    )
    featurecounts_changes: Dict[str, Any] = Field(
        default={},
        description="FeatureCounts参数修改：-s, -p, -M, -O, -Q等FeatureCounts特有参数。当用户明确提到FeatureCounts或这些参数时必须使用此字段！"
    )
    modification_reason: str = Field(
        default="", 
        description="修改原因说明"
    )
    validation_notes: List[str] = Field(
        default=[], 
        description="参数验证提示"
    )


async def modify_node(state: AgentState) -> Dict[str, Any]:
    """
    修改节点 - 处理用户的配置修改请求
    
    工作流程：
    1. 接收用户修改需求（自然语言）
    2. 使用LLM解析为结构化修改
    3. 验证并应用修改到当前状态
    4. 更新所有相关配置字段
    5. 直接返回到User Confirm节点展示更新后的配置
    
    设计理念：
    - Modify Node 直接修改状态，无需重新生成
    - Prepare Node 只在初始配置时使用一次
    - 所有后续修改都通过 Modify Node 完成
    """
    logger = get_logger("rna.nodes.modify")
    logger.info("配置修改节点启动")
    
    # 获取修改需求 - 在state.input获取用户输入
    raw_input = state.input or ""
    
    # 获取当前配置
    current_nextflow = state.nextflow_config or {}
    current_resource = state.resource_config or {}
    current_fastp = state.fastp_params or {}
    current_star = state.star_params or {}
    current_hisat2 = state.hisat2_params or {}
    current_featurecounts = state.featurecounts_params or {}
    
    logger.info(f"用户修改需求: {raw_input}")
    logger.debug(
        "当前配置概要 | nextflow=%d 资源进程=%d fastp=%d star=%d hisat2=%d featurecounts=%d",
        len(current_nextflow), len(current_resource), len(current_fastp), len(current_star), len(current_hisat2), len(current_featurecounts)
    )
    
    # 分析当前执行上下文 - 获取工具选择
    current_step = state.current_step or ""
    completed_steps = state.completed_steps or []
    nextflow_config = state.nextflow_config or {}

    # 获取所有工具配置
    qc_tool = nextflow_config.get("qc_tool", "fastp")  # 默认FastP
    align_tool = nextflow_config.get("align_tool", "star")  # 默认STAR
    quant_tool = nextflow_config.get("quant_tool", "featurecounts")  # 默认FeatureCounts

    logger.info(f"执行上下文 | 当前步骤={current_step} 已完成={completed_steps}")
    logger.info(f"工具配置 | 质控={qc_tool} 比对={align_tool} 定量={quant_tool}")
    
    # 构建上下文感知的LLM提示
    context_analysis = []
    if current_step:
        context_analysis.append(f"- 当前正在执行: {current_step}")
    if completed_steps:
        context_analysis.append(f"- 已完成步骤: {', '.join(completed_steps)}")
    context_analysis.append(f"- 配置的质控工具: {qc_tool}")
    context_analysis.append(f"- 配置的比对工具: {align_tool}")
    context_analysis.append(f"- 配置的定量工具: {quant_tool}")
    
    context_info = "\n".join(context_analysis) if context_analysis else "- 尚未开始执行"

    user_prompt = f"""🎯 **当前执行上下文**：
{context_info}

用户修改需求：
{raw_input}

请解析修改需求并返回结构化的修改内容。
"""

    # 调用LLM解析修改需求
    try:
        llm = get_shared_llm()
        llm_with_structure = llm.with_structured_output(ModifyRequest)
        
        # 构建LangGraph标准消息格式
        messages = [
            {"role": "system", "content": MODIFY_NODE_PROMPT},
            {"role": "user", "content": user_prompt}
        ]
        
        logger.info("解析修改需求（调用LLM）...")
        modify_request = await llm_with_structure.ainvoke(messages)
        try:
            log_llm_preview(logger, "modify", modify_request)
        except Exception:
            pass
        
        # 应用修改
        updated_nextflow = current_nextflow.copy()
        updated_resource = current_resource.copy()
        updated_fastp = current_fastp.copy()
        updated_star = current_star.copy()
        updated_hisat2 = current_hisat2.copy()
        updated_featurecounts = current_featurecounts.copy()
        
        # 应用Nextflow配置修改
        if modify_request.nextflow_changes:
            logger.info("应用Nextflow配置修改")
            for key, value in modify_request.nextflow_changes.items():
                old_value = updated_nextflow.get(key, "未设置")
                updated_nextflow[key] = value
                logger.debug(f"nextflow.{key}: {old_value} → {value}")
        
        # 应用资源配置修改
        if modify_request.resource_changes:
            logger.info("应用资源配置修改")
            for process, changes in modify_request.resource_changes.items():
                if process not in updated_resource:
                    updated_resource[process] = {}
                for key, value in changes.items():
                    old_value = updated_resource[process].get(key, "未设置")
                    updated_resource[process][key] = value
                    logger.debug(f"{process}.{key}: {old_value} → {value}")
        
        # 应用FastP参数修改（统一键名策略：仅接受精确键名，忽略未知键）
        if modify_request.fastp_changes:
            logger.info("应用FastP参数修改")

            allowed_keys = {
                "qualified_quality_phred", "unqualified_percent_limit", "n_base_limit", "length_required",
                "adapter_trimming", "quality_filtering", "length_filtering",
                "phred64", "reads_to_process", "fix_mgi_id", "detect_adapter_for_pe",
                "trim_front1", "trim_tail1", "max_len1", "trim_front2", "trim_tail2", "max_len2",
                "trim_poly_g", "poly_g_min_len", "disable_trim_poly_g", "trim_poly_x", "poly_x_min_len",
                "cut_front", "cut_tail", "cut_right", "cut_window_size", "cut_mean_quality",
                "cut_front_window_size", "cut_front_mean_quality", "cut_tail_window_size", "cut_tail_mean_quality", "cut_right_window_size", "cut_right_mean_quality",
                "average_qual", "disable_length_filtering", "length_limit", "low_complexity_filter", "complexity_threshold",
                "correction", "overlap_len_require", "overlap_diff_limit", "overlap_diff_percent_limit",
                "overrepresentation_analysis", "overrepresentation_sampling"
            }

            def _to_bool(v):
                if isinstance(v, bool):
                    return v
                if isinstance(v, str):
                    return v.strip().lower() in {"1", "true", "yes", "y", "on"}
                if isinstance(v, (int, float)):
                    return bool(v)
                return False

            for key, value in modify_request.fastp_changes.items():
                if key not in allowed_keys:
                    logger.warning(f"跳过未知FastP键: {key}")
                    continue
                if key in {
                    "adapter_trimming", "quality_filtering", "length_filtering",
                    "phred64", "fix_mgi_id", "trim_poly_g", "disable_trim_poly_g",
                    "trim_poly_x", "cut_front", "cut_tail", "cut_right",
                    "low_complexity_filter", "correction", "overrepresentation_analysis",
                    "detect_adapter_for_pe",
                }:
                    value = _to_bool(value)
                old_value = updated_fastp.get(key, "未设置")
                updated_fastp[key] = value
                logger.debug(f"fastp.{key}: {old_value} → {value}")
        
        # 应用STAR参数修改
        if modify_request.star_changes:
            logger.info("应用STAR参数修改")
            
            star_allowed_keys = {
                "outSAMtype", "outSAMunmapped", "outSAMattributes",
                "outFilterMultimapNmax", "alignSJoverhangMin", "alignSJDBoverhangMin", 
                "outFilterMismatchNmax", "outFilterMismatchNoverReadLmax",
                "alignIntronMin", "alignIntronMax", "alignMatesGapMax", "quantMode", "twopassMode",
                "limitBAMsortRAM", "outBAMsortingThreadN", "genomeLoad", "outFileNamePrefix",
                "readFilesCommand", "outReadsUnmapped", "outFilterIntronMotifs", 
                "outSAMstrandField", "outFilterType", "sjdbGTFfile", "sjdbOverhang", 
                "chimSegmentMin", "chimOutType", "chimMainSegmentMultNmax"
            }
            
            for key, value in modify_request.star_changes.items():
                if key not in star_allowed_keys:
                    logger.warning(f"跳过未知STAR键: {key}")
                    continue
                old_value = updated_star.get(key, "未设置")
                updated_star[key] = value
                logger.debug(f"star.{key}: {old_value} → {value}")
        
        # 应用HISAT2参数修改
        if modify_request.hisat2_changes:
            logger.info("应用HISAT2参数修改")
            
            hisat2_allowed_keys = {
                "--mp", "--rdg", "--rfg", "--score-min", "--ma", "--np", "--sp", "--no-mixed", "--no-discordant",
                "--gbar", "--ignore-quals", "--nofw", "--norc", "--end-to-end", "--local", "--very-fast",
                "--fast", "--sensitive", "--very-sensitive", "--very-fast-local", "--fast-local", 
                "--sensitive-local", "--very-sensitive-local", "-N", "-L", "-i", "--n-ceil",
                "-D", "-R", "--dpad", "--gbar", "--ignore-quals", "--nofw", "--norc", "--no-1mm-upfront",
                "-k", "-a", "--time", "--un", "--al", "--un-conc", "--al-conc", "--summary-file",
                "--new-summary", "--quiet", "--met-file", "--met-stderr", "--met", "--no-head",
                "--no-sq", "--rg-id", "--rg", "--omit-sec-seq", "--sam-no-qname-trunc", "--xeq",
                "--soft-clipped-unmapped-tlen", "--sam-append-comment", "--reorder", "--mm",
                "--qc-filter", "--seed", "--non-deterministic", "--remove-chrname-prefix", "--add-chrname-prefix"
            }
            
            for key, value in modify_request.hisat2_changes.items():
                if key not in hisat2_allowed_keys:
                    logger.warning(f"跳过未知HISAT2键: {key}")
                    continue
                # 处理布尔类型参数
                if key in {"--no-mixed", "--no-discordant", "--ignore-quals", "--nofw", "--norc", 
                          "--end-to-end", "--local", "--no-1mm-upfront", "--time", "--new-summary", 
                          "--quiet", "--met-stderr", "--no-head", "--no-sq", "--omit-sec-seq", 
                          "--sam-no-qname-trunc", "--xeq", "--soft-clipped-unmapped-tlen", 
                          "--sam-append-comment", "--reorder", "--qc-filter", "--non-deterministic", 
                          "--remove-chrname-prefix", "--add-chrname-prefix"}:
                    if isinstance(value, str):
                        value = value.strip().lower() in {"1", "true", "yes", "y", "on"}
                    elif isinstance(value, (int, float)):
                        value = bool(value)
                old_value = updated_hisat2.get(key, "未设置")
                updated_hisat2[key] = value
                logger.debug(f"hisat2.{key}: {old_value} → {value}")
        
        # 应用FeatureCounts参数修改
        if modify_request.featurecounts_changes:
            logger.info("应用FeatureCounts参数修改")
            
            fc_allowed_keys = {
                "-s", "-p", "-B", "-C", "-t", "-g", "-M", "-O", "--fraction", "-Q",
                "--minOverlap", "--fracOverlap", "-f", "-J",
                "-a", "-F", "--primary", "--ignoreDup", "--splitOnly", "--nonSplitOnly", 
                "--largestOverlap", "--readShiftType", "--readShiftSize", "-R", 
                "--readExtension5", "--readExtension3", "--read2pos", "--countReadPairs",
                "--donotsort", "--byReadGroup", "--extraAttributes"
            }
            
            for key, value in modify_request.featurecounts_changes.items():
                if key not in fc_allowed_keys:
                    logger.warning(f"跳过未知FeatureCounts键: {key}")
                    continue
                # 处理布尔类型参数
                if key in {"-p", "-B", "-C", "-M", "-O", "--fraction", "-f", "-J", 
                          "--primary", "--ignoreDup", "--splitOnly", "--nonSplitOnly", "--largestOverlap"}:
                    value = _to_bool(value)
                old_value = updated_featurecounts.get(key, "未设置")
                updated_featurecounts[key] = value
                logger.debug(f"featurecounts.{key}: {old_value} → {value}")
        
        # 显示验证提示
        if modify_request.validation_notes:
            for note in modify_request.validation_notes:
                logger.warning(f"参数验证提示: {note}")
        
        # 记录修改历史
        modification_history = getattr(state, 'modification_history', []) or []
        modification_record = {
            "timestamp": datetime.now().isoformat(),
            "raw_input": raw_input,
            "changes": {
                "nextflow": modify_request.nextflow_changes,
                "resource": modify_request.resource_changes,
                "fastp": modify_request.fastp_changes,
                "star": modify_request.star_changes,
                "featurecounts": modify_request.featurecounts_changes
            },
            "reason": modify_request.modification_reason
        }
        modification_record["changes"]["hisat2"] = modify_request.hisat2_changes
        modification_history.append(modification_record)
        
        logger.info("配置修改完成")
        logger.info(f"修改原因: {modify_request.modification_reason}")

        # 直接更新 Prepare 节点生成的参数文件（不做版本管理）
        try:
            # 目标目录：统一使用 state.results_dir
            target_results_dir = getattr(state, "results_dir", "")
            if target_results_dir:
                params_dir = Path(target_results_dir) / "params"
                params_dir.mkdir(parents=True, exist_ok=True)

                # 查找最新的 prepare_params_*.json；若不存在则使用固定文件名
                prepare_files = sorted(params_dir.glob("prepare_params_*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
                if prepare_files:
                    prepare_path = prepare_files[0]
                else:
                    prepare_path = params_dir / "prepare_params.json"

                # 载入原文件（如存在），仅更新 nextflow_config 与 resource_config
                payload: Dict[str, Any] = {}
                if prepare_path.exists():
                    try:
                        with open(prepare_path, "r", encoding="utf-8") as f:
                            payload = json.load(f) or {}
                    except Exception:
                        payload = {}

                payload["nextflow_config"] = updated_nextflow
                payload["resource_config"] = updated_resource

                with open(prepare_path, "w", encoding="utf-8") as f:
                    json.dump(payload, f, indent=2, ensure_ascii=False)

                logger.info(f"已更新 Prepare 参数文件: {prepare_path}")
            else:
                logger.warning("无法定位 results_dir，跳过更新 Prepare 参数文件")
        except Exception as e:
            logger.warning(f"更新 Prepare 参数文件失败（不影响流程）: {e}")

        logger.info("返回到确认节点查看更新后的配置")
        
        # 返回更新后的状态
        return {
            "success": True,
            # 更新配置
            "nextflow_config": updated_nextflow,
            "resource_config": updated_resource,
            "fastp_params": updated_fastp,
            "star_params": updated_star,
            "hisat2_params": updated_hisat2,
            "featurecounts_params": updated_featurecounts,
            
            # 记录修改处理结果
            "modify_results": {
                "original_input": raw_input,
                "parsed_changes": {
                    "nextflow_config": modify_request.nextflow_changes,
                    "resource_config": modify_request.resource_changes,
                    "fastp_params": modify_request.fastp_changes,
                    "star_params": modify_request.star_changes,
                    "hisat2_params": modify_request.hisat2_changes,
                    "featurecounts_params": modify_request.featurecounts_changes
                },
                "applied": True
            },
            
            # 保存修改历史
            "modification_history": modification_history,
            
            # 配置理由保持原样，避免重复插入修改说明
            "config_reasoning": state.config_reasoning,
            
            # 状态和响应
            "response": f"✅ 已应用配置修改：{modify_request.modification_reason}",
            "status": "success"  # 修改成功状态
        }
        
    except Exception as e:
        logger.error(f"修改解析失败: {str(e)}", exc_info=True)
        logger.info("返回到用户确认节点")
        
        return {
            "success": False,
            "response": f"❌ 修改解析失败: {str(e)}",
            "status": "failed",
            "modify_results": {}  # 清空修改结果
        }
