"""
Modify Node - 智能配置修改节点
负责解析用户修改需求并更新所有相关配置参数
"""
from typing import Dict, Any, List
from datetime import datetime
from pydantic import BaseModel, Field
from ..state import AgentState
from ..core import get_shared_llm
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
    
    print(f"\n{'='*60}")
    print(f"🔧 **配置修改节点**")
    print(f"{'='*60}")
    
    # 获取修改需求
    modify_requirements = state.modify_requirements or {}
    raw_input = modify_requirements.get("raw_input", "")
    
    # 获取当前配置
    current_nextflow = state.nextflow_config or {}
    current_resource = state.resource_config or {}
    current_fastp = state.fastp_params or {}
    current_star = state.star_params or {}
    current_featurecounts = state.featurecounts_params or {}
    
    print(f"\n📝 用户修改需求: {raw_input}")
    print(f"\n📋 当前配置概览:")
    print(f"   - Nextflow配置: {len(current_nextflow)} 项")
    print(f"   - 资源配置: {len(current_resource)} 个进程")
    print(f"   - FastP参数: {len(current_fastp)} 项")
    print(f"   - STAR参数: {len(current_star)} 项")
    print(f"   - FeatureCounts参数: {len(current_featurecounts)} 项")
    
    # 构建LLM提示
    system_prompt = """你是RNA-seq分析配置专家。请解析用户的修改需求，将其转换为具体的参数修改。

‼️ **必须遵守的字段选择规则**：
1. 如果用户提到"STAR"、"outFilterMultimapNmax"、"twopassMode"等STAR相关参数 → 使用star_changes字段
2. 如果用户提到"FeatureCounts"、"-s"、"-p"、"-M"、"-Q"等FeatureCounts相关参数 → 使用featurecounts_changes字段
3. 如果用户提到"FastP"、"qualified_quality_phred"、"length_required"等FastP相关参数 → 使用fastp_changes字段

‼️ **绝对禁止**：不要说参数"不在配置范围内"！用户当前提供了完整的STAR和FeatureCounts参数，你必须使用对应的字段！

严格要求：请使用下方【精确键名】返回修改，禁止使用任何别名或同义词；布尔值请使用 true/false，数值使用数字。

【Nextflow配置参数（键名必须精确）】
- species, genome_version, qc_tool, align_tool, quant_tool, paired_end,
- run_download_genome, run_build_star_index, run_build_hisat2_index

【资源配置参数（按进程）】
- 每个进程键名与字段：{"<process>": {"cpus": <int>, "memory": "<GB字符串>"}}
- 进程：prepare_star_index, prepare_hisat2_index, run_alignment, run_quality_control, run_quantification, download_genome_fasta, download_genome_gtf

【FastP参数（键名必须精确）】
- qualified_quality_phred, unqualified_percent_limit, n_base_limit, length_required,
- adapter_trimming, quality_filtering, length_filtering,
- phred64, reads_to_process, fix_mgi_id, detect_adapter_for_pe,
- trim_front1, trim_tail1, max_len1, trim_front2, trim_tail2, max_len2,
- trim_poly_g, poly_g_min_len, disable_trim_poly_g, trim_poly_x, poly_x_min_len,
- cut_front, cut_tail, cut_right, cut_window_size, cut_mean_quality,
- cut_front_window_size, cut_front_mean_quality, cut_tail_window_size, cut_tail_mean_quality, cut_right_window_size, cut_right_mean_quality,
- average_qual, disable_length_filtering, length_limit, low_complexity_filter, complexity_threshold,
- correction, overlap_len_require, overlap_diff_limit, overlap_diff_percent_limit,
- overrepresentation_analysis, overrepresentation_sampling

【STAR参数（键名必须精确）】
- outSAMtype, outSAMunmapped, outSAMattributes,
- outFilterMultimapNmax, alignSJoverhangMin, alignSJDBoverhangMin, outFilterMismatchNmax, outFilterMismatchNoverReadLmax,
- alignIntronMin, alignIntronMax, alignMatesGapMax, quantMode, twopassMode,
- limitBAMsortRAM, outBAMsortingThreadN, genomeLoad, outFileNamePrefix,
- readFilesCommand, outReadsUnmapped, outFilterIntronMotifs, outSAMstrandField,
- outFilterType, sjdbGTFfile, sjdbOverhang, chimSegmentMin, chimOutType, chimMainSegmentMultNmax

【FeatureCounts参数（键名必须精确）】
- -s, -p, -B, -C, -t, -g, -M, -O, --fraction, -Q,
- --minOverlap, --fracOverlap, -f, -J,
- -a, -F, --primary, --ignoreDup, --splitOnly, --nonSplitOnly, --largestOverlap,
- --readShiftType, --readShiftSize, -R, --readExtension5, --readExtension3,
- --read2pos, --countReadPairs, --donotsort, --byReadGroup, --extraAttributes

⚠️ **关键参数选择规则**：
1. **质量相关参数** → 使用 fastp_changes：如"质量阈值"、"qualified_quality_phred"、"length_required"
2. **比对相关参数** → 使用 star_changes：如"多重比对"、"两遍模式"、"outFilterMultimapNmax"、"twopassMode"  
3. **计数相关参数** → 使用 featurecounts_changes：如"链特异性"、"双端模式"、"-s"、"-p"、"-M"
4. **线程/CPU资源** → 使用 resource_changes：如"线程数"、"CPU核心"、"runThreadN"、"-T"参数
5. **流程配置** → 使用 nextflow_changes：物种、基因组版本、工具选择

⚠️ **重要提醒**：用户明确提到具体工具参数时，必须使用对应的工具参数字段！

请分析用户需求，优先使用工具专用参数字段，返回需要修改的参数。只修改用户明确要求的部分，保持其他配置不变，并严格使用上述精确键名。
"""

    user_prompt = f"""当前配置状态：

Nextflow配置：
{json.dumps(current_nextflow, indent=2, ensure_ascii=False)}

资源配置：
{json.dumps(current_resource, indent=2, ensure_ascii=False)}

FastP参数：
{json.dumps(current_fastp, indent=2, ensure_ascii=False)}

STAR参数：
{json.dumps(current_star, indent=2, ensure_ascii=False)}

FeatureCounts参数：
{json.dumps(current_featurecounts, indent=2, ensure_ascii=False)}

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
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ]
        
        print(f"\n🤖 正在解析修改需求...")
        modify_request = llm_with_structure.invoke(messages)
        
        # 应用修改
        updated_nextflow = current_nextflow.copy()
        updated_resource = current_resource.copy()
        updated_fastp = current_fastp.copy()
        updated_star = current_star.copy()
        updated_featurecounts = current_featurecounts.copy()
        
        # 应用Nextflow配置修改
        if modify_request.nextflow_changes:
            print(f"\n📦 应用Nextflow配置修改:")
            for key, value in modify_request.nextflow_changes.items():
                old_value = updated_nextflow.get(key, "未设置")
                updated_nextflow[key] = value
                print(f"   - {key}: {old_value} → {value}")
        
        # 应用资源配置修改
        if modify_request.resource_changes:
            print(f"\n💻 应用资源配置修改:")
            for process, changes in modify_request.resource_changes.items():
                if process not in updated_resource:
                    updated_resource[process] = {}
                for key, value in changes.items():
                    old_value = updated_resource[process].get(key, "未设置")
                    updated_resource[process][key] = value
                    print(f"   - {process}.{key}: {old_value} → {value}")
        
        # 应用FastP参数修改（统一键名策略：仅接受精确键名，忽略未知键）
        if modify_request.fastp_changes:
            print(f"\n🧬 应用FastP参数修改:")

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
                    print(f"   - 跳过未知键: {key}")
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
                print(f"   - {key}: {old_value} → {value}")
        
        # 应用STAR参数修改
        if modify_request.star_changes:
            print(f"\n⭐ 应用STAR参数修改:")
            
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
                    print(f"   - 跳过未知STAR键: {key}")
                    continue
                old_value = updated_star.get(key, "未设置")
                updated_star[key] = value
                print(f"   - {key}: {old_value} → {value}")
        
        # 应用FeatureCounts参数修改
        if modify_request.featurecounts_changes:
            print(f"\n📊 应用FeatureCounts参数修改:")
            
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
                    print(f"   - 跳过未知FeatureCounts键: {key}")
                    continue
                # 处理布尔类型参数
                if key in {"-p", "-B", "-C", "-M", "-O", "--fraction", "-f", "-J", 
                          "--primary", "--ignoreDup", "--splitOnly", "--nonSplitOnly", "--largestOverlap"}:
                    value = _to_bool(value)
                old_value = updated_featurecounts.get(key, "未设置")
                updated_featurecounts[key] = value
                print(f"   - {key}: {old_value} → {value}")
        
        # 显示验证提示
        if modify_request.validation_notes:
            print(f"\n⚠️ 参数验证提示:")
            for note in modify_request.validation_notes:
                print(f"   - {note}")
        
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
        modification_history.append(modification_record)
        
        print(f"\n✅ 配置修改完成！")
        print(f"💭 修改原因: {modify_request.modification_reason}")
        print(f"\n🔄 返回到确认节点查看更新后的配置...")
        
        # 返回更新后的状态
        return {
            "success": True,
            # 更新配置
            "nextflow_config": updated_nextflow,
            "resource_config": updated_resource,
            "fastp_params": updated_fastp,
            "star_params": updated_star,
            "featurecounts_params": updated_featurecounts,
            
            # 更新修改需求（记录已应用）
            "modify_requirements": {
                "raw_input": raw_input,
                "parsed_changes": {
                    "nextflow_config": modify_request.nextflow_changes,
                    "resource_config": modify_request.resource_changes,
                    "fastp_params": modify_request.fastp_changes,
                    "star_params": modify_request.star_changes,
                    "featurecounts_params": modify_request.featurecounts_changes
                },
                "applied": True
            },
            
            # 保存修改历史
            "modification_history": modification_history,
            
            # 配置理由保持原样（不再在理由中插入修改说明，避免重复展示）
            "config_reasoning": state.config_reasoning,
            
            # 状态和响应
            "response": f"✅ 已应用配置修改：{modify_request.modification_reason}",
            "status": "success"  # 修改成功状态
        }
        
    except Exception as e:
        print(f"\n❌ 修改解析失败: {str(e)}")
        print(f"🔄 返回到用户确认节点...")
        
        return {
            "success": False,
            "response": f"❌ 修改解析失败: {str(e)}",
            "status": "failed",
            "modify_requirements": {}  # 清空修改需求
        }
