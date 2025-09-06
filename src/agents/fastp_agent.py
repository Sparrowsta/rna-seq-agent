"""
FastP Agent - RNA-seq质量控制和预处理Agent
支持智能配置、结果分析和参数优化
"""

import json
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field

from langchain_core.tools import tool
from langchain_core.messages import HumanMessage

# 导入项目配置
from ..config import get_tools_config
from ..core import LLMManager


@dataclass
class FastpConfig:
    """Fastp配置参数类"""
    
    # 基本参数
    input_files: List[str]  # 输入FASTQ文件列表
    output_dir: str  # 输出目录
    sample_name: str  # 样本名称
    
    # fastp参数
    adapter_trimming: bool = True  # 自动adapter trimming
    quality_filtering: bool = True  # 质量过滤
    length_filtering: bool = True  # 长度过滤
    
    # 质量参数
    qualified_quality_phred: int = 20  # 质量阈值
    unqualified_percent_limit: int = 40  # 低质量base比例限制
    n_base_limit: int = 5  # N碱基数量限制
    
    # 长度参数
    length_required: int = 15  # 最短长度要求
    
    # 输出参数
    html_report: bool = True  # 生成HTML报告
    json_report: bool = True  # 生成JSON报告
    
    fastp_cpus: Optional[int] = None  # 若通过 Nextflow 运行，作为 process.cpus 传入

    # 额外 fastp 高级参数（按需传递）
    # 输入质量编码与读取控制
    phred64: bool = False                 # -6/--phred64 输入为phred64质量
    reads_to_process: Optional[int] = None  # --reads_to_process 仅处理前N条reads
    fix_mgi_id: bool = False              # --fix_mgi_id 修复MGI测序ID

    # PE adapter 自动检测
    detect_adapter_for_pe: Optional[bool] = None  # 若为None按自动逻辑；显式True/False可覆盖

    # 前后端定长修剪与最大长度
    trim_front1: Optional[int] = None     # -f/--trim_front1
    trim_tail1: Optional[int] = None      # -t/--trim_tail1
    max_len1: Optional[int] = None        # -b/--max_len1
    trim_front2: Optional[int] = None     # -F/--trim_front2
    trim_tail2: Optional[int] = None      # -T/--trim_tail2
    max_len2: Optional[int] = None        # -B/--max_len2

    # polyG / polyX 修剪
    trim_poly_g: Optional[bool] = None    # -g/--trim_poly_g（Illumina NextSeq/NovaSeq常见）
    poly_g_min_len: Optional[int] = None  # --poly_g_min_len
    disable_trim_poly_g: Optional[bool] = None  # -G/--disable_trim_poly_g
    trim_poly_x: Optional[bool] = None    # -x/--trim_poly_x
    poly_x_min_len: Optional[int] = None  # --poly_x_min_len

    # 滑窗切除与门限
    cut_front: Optional[int] = None       # -5/--cut_front
    cut_tail: Optional[int] = None        # -3/--cut_tail
    cut_right: Optional[int] = None       # -r/--cut_right
    cut_window_size: Optional[int] = None # -W/--cut_window_size
    cut_mean_quality: Optional[int] = None# -M/--cut_mean_quality
    cut_front_window_size: Optional[int] = None
    cut_front_mean_quality: Optional[int] = None
    cut_tail_window_size: Optional[int] = None
    cut_tail_mean_quality: Optional[int] = None
    cut_right_window_size: Optional[int] = None
    cut_right_mean_quality: Optional[int] = None

    # 质量/长度过滤细化
    average_qual: Optional[int] = None    # -e/--average_qual
    disable_length_filtering: Optional[bool] = None  # -L/--disable_length_filtering（与 length_filtering 相反）
    length_limit: Optional[int] = None    # --length_limit
    low_complexity_filter: Optional[bool] = None  # -y/--low_complexity_filter
    complexity_threshold: Optional[int] = None    # -Y/--complexity_threshold

    # PE 重叠校正与检测
    correction: Optional[bool] = None     # -c/--correction（仅PE）
    overlap_len_require: Optional[int] = None
    overlap_diff_limit: Optional[int] = None
    overlap_diff_percent_limit: Optional[int] = None

    # 过表达序列分析
    overrepresentation_sampling: Optional[int] = None  # -P/--overrepresentation_sampling


@dataclass
class FastpResult:
    """Fastp执行结果类"""
    
    # 执行状态
    success: bool
    exit_code: int
    error_message: str = ""
    
    # 输出文件
    output_files: List[str] = field(default_factory=list)
    html_report: str = ""
    json_report: str = ""
    
    # 基本统计信息
    total_reads_before: int = 0
    total_reads_after: int = 0
    total_bases_before: int = 0
    total_bases_after: int = 0
    q20_bases_before: int = 0
    q20_bases_after: int = 0
    q30_bases_before: int = 0
    q30_bases_after: int = 0
    
    # 测序信息
    sequencing_type: str = ""  # single end或paired end
    read_length_before: float = 0  # 平均长度
    read_length_after: float = 0
    
    # GC含量
    gc_content_before: float = 0.0
    gc_content_after: float = 0.0
    
    # 过滤统计
    passed_filter_reads: int = 0
    low_quality_reads: int = 0
    too_many_n_reads: int = 0
    too_short_reads: int = 0
    too_long_reads: int = 0
    
    # 接头信息
    adapter_trimmed_reads: int = 0
    adapter_trimmed_bases: int = 0
    detected_adapters: List[str] = field(default_factory=list)
    
    # 重复率
    duplication_rate: float = 0.0
    
    # 过度表达序列
    overrepresented_sequences: int = 0
    top_overrepresented: List[Dict[str, Any]] = field(default_factory=list)
    
    # 质量评估
    mean_quality_before: float = 0.0
    mean_quality_after: float = 0.0
    
    def __post_init__(self):
        # 不再需要手动初始化，field(default_factory=list)已处理
        pass


class FastpAgent:
    """
    FastP批次质量控制与参数优化Agent
    
    统一的批次处理架构，支持：
    - 批次样本处理（单样本作为批次大小1的特例）
    - 智能参数优化建议（基于批次结果分析）
    - 容器化Nextflow执行环境
    - LLM智能质量评估和优化策略
    """
    
    def __init__(self):
        """初始化FastP批次处理Agent"""
        self.config_manager = get_tools_config()
        # 初始化LLM管理器
        self.llm_manager = LLMManager(self.config_manager.settings)
        
    def validate_config(self, config: FastpConfig) -> Tuple[bool, List[str]]:
        """
        验证配置参数的有效性
        
        Args:
            config: FastP配置对象
            
        Returns:
            (is_valid, error_messages): 验证结果和错误信息列表
        """
        errors = []
        
        # 验证输入文件
        if not config.input_files:
            errors.append("输入文件列表不能为空")
        else:
            for file_path in config.input_files:
                file_obj = Path(file_path)
                if not file_obj.exists():
                    errors.append(f"输入文件不存在: {file_path}")
                elif file_obj.stat().st_size == 0:
                    errors.append(f"输入文件为空: {file_path}")
        
        # 验证输出目录
        if not config.output_dir:
            errors.append("输出目录不能为空")
        else:
            output_path = Path(config.output_dir)
            try:
                output_path.mkdir(parents=True, exist_ok=True)
            except PermissionError:
                errors.append(f"无权限创建输出目录: {config.output_dir}")
        
        # 验证样本名称
        if not config.sample_name:
            errors.append("样本名称不能为空")
        elif not config.sample_name.replace('_', '').replace('-', '').isalnum():
            errors.append("样本名称只能包含字母、数字、下划线和连字符")
        
        # 验证质量参数
        if config.qualified_quality_phred < 0 or config.qualified_quality_phred > 40:
            errors.append("质量阈值必须在0-40之间")
        
        if config.unqualified_percent_limit < 0 or config.unqualified_percent_limit > 100:
            errors.append("低质量base比例限制必须在0-100之间")
        
        if config.length_required < 1:
            errors.append("最短长度要求必须大于0")
        
        return len(errors) == 0, errors
    
    def _parse_json_report_with_path(self, json_path: str, result: FastpResult) -> FastpResult:
        """使用路径解析fastp JSON报告，提取详细统计信息"""
        try:
            import json
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # 提取summary统计信息
            summary = data.get('summary', {})
            
            # 测序类型
            result.sequencing_type = summary.get('sequencing', 'unknown')
            
            # 处理前统计
            before_filtering = summary.get('before_filtering', {})
            result.total_reads_before = before_filtering.get('total_reads', 0)
            result.total_bases_before = before_filtering.get('total_bases', 0)
            result.q20_bases_before = before_filtering.get('q20_bases', 0)
            result.q30_bases_before = before_filtering.get('q30_bases', 0)
            result.read_length_before = before_filtering.get('read1_mean_length', 0)
            result.gc_content_before = before_filtering.get('gc_content', 0.0)
            
            # 处理后统计
            after_filtering = summary.get('after_filtering', {})
            result.total_reads_after = after_filtering.get('total_reads', 0)
            result.total_bases_after = after_filtering.get('total_bases', 0)
            result.q20_bases_after = after_filtering.get('q20_bases', 0)
            result.q30_bases_after = after_filtering.get('q30_bases', 0)
            result.read_length_after = after_filtering.get('read1_mean_length', 0)
            result.gc_content_after = after_filtering.get('gc_content', 0.0)
            
            # 过滤统计
            filtering_result = data.get('filtering_result', {})
            result.passed_filter_reads = filtering_result.get('passed_filter_reads', 0)
            result.low_quality_reads = filtering_result.get('low_quality_reads', 0)
            result.too_many_n_reads = filtering_result.get('too_many_N_reads', 0)
            result.too_short_reads = filtering_result.get('too_short_reads', 0)
            result.too_long_reads = filtering_result.get('too_long_reads', 0)
            
            # 重复率
            duplication = data.get('duplication', {})
            result.duplication_rate = duplication.get('rate', 0.0)
            
            # 接头切除信息
            adapter_cutting = data.get('adapter_cutting', {})
            result.adapter_trimmed_reads = adapter_cutting.get('adapter_trimmed_reads', 0)
            result.adapter_trimmed_bases = adapter_cutting.get('adapter_trimmed_bases', 0)
            
            # 检测到的接头序列
            if 'read1_adapter_sequence' in adapter_cutting:
                result.detected_adapters.append(adapter_cutting['read1_adapter_sequence'])
            if 'read2_adapter_sequence' in adapter_cutting:
                result.detected_adapters.append(adapter_cutting['read2_adapter_sequence'])
            
            # 过度表达序列统计
            overrep = data.get('read1_after_filtering', {}).get('overrepresented_sequences', {})
            if not overrep:  # 如果read1没有，尝试summary级别
                overrep = data.get('overrepresented_sequences', {})
            
            if overrep:
                result.overrepresented_sequences = len(overrep)
                # 提取前5个过度表达序列
                top_seqs = []
                for seq, count in sorted(overrep.items(), key=lambda x: x[1] if isinstance(x[1], int) else 0, reverse=True)[:5]:
                    top_seqs.append({'sequence': seq, 'count': count})
                result.top_overrepresented = top_seqs
                
        except Exception as e:
            # JSON解析失败不影响主要功能，记录错误即可
            result.error_message = f"解析JSON报告失败: {str(e)}"
        
        return result
    
    def _parse_llm_json_response(self, response_text: str) -> Dict[str, Any]:
        """解析LLM返回的JSON响应，处理代码块标记"""
        # 移除代码块标记
        if response_text.startswith("```"):
            start_idx = response_text.find("{")
            end_idx = response_text.rfind("}") + 1
            if start_idx != -1 and end_idx > start_idx:
                response_text = response_text[start_idx:end_idx]
        
        import json
        return json.loads(response_text)
    
    def _get_llm_group_summary_and_params(
        self,
        samples: List[Dict[str, Any]],
        base_params: Dict[str, Any],
        history: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """
        统一的LLM分析方法：将所有样本结果提交给LLM，获取总结与参数建议
        
        返回结构：
        {
          "suggested_params": {k: v, ...},   # 仅需变更的键；作为覆盖项
          "reasoning": "...",               # 参数推荐理由（中文）
          "group_summary": "...",           # 面向用户的整体质量总结（中文）
          "sample_assessments": {"sample_id": "质量评级", ...}  # 各样本质量评级
        }
        """
        # 提取完整的样本统计信息供LLM分析
        compact = []
        for s in samples:
            if not s.get("success"):
                compact.append({
                    "sample_id": s.get("sample_id"),
                    "success": False,
                    "error": s.get("error", "未知错误")
                })
                continue
                
            stats = s.get("statistics", {})
            
            # 从results中已有的统计数据构建完整信息
            sample_info = {
                "sample_id": s.get("sample_id"),
                "success": True,
                # 基本统计 - 这些在 run_batch 中已经计算好了
                "reads_before": stats.get("reads_before", 0),
                "reads_after": stats.get("reads_after", 0),
                "reads_retention": stats.get("reads_retention", 0.0),
                "q30_rate": stats.get("q30_rate", 0.0),
                # 从 FastpResult 中提取的额外信息
                "gc_content_before": stats.get("gc_content_before", 0.0),
                "gc_content_after": stats.get("gc_content_after", 0.0),
                "mean_length_before": stats.get("mean_length_before", 0.0),
                "mean_length_after": stats.get("mean_length_after", 0.0),
                "duplication_rate": stats.get("duplication_rate", 0.0),
                "adapter_trimmed_reads": stats.get("adapter_trimmed_reads", 0),
                "adapter_trimmed_percent": stats.get("adapter_trimmed_percent", 0.0),
                # 过滤原因统计
                "low_quality_reads": stats.get("low_quality_reads", 0),
                "too_short_reads": stats.get("too_short_reads", 0),
                "too_many_n_reads": stats.get("too_many_n_reads", 0),
            }
            
            compact.append(sample_info)

        import json as _json
        sys_prompt = (
            "你是资深生信质控工程师。请基于一组 fastp 结果给出整体总结并推荐参数。"
            " 输出严格 JSON（不含代码块、注释），键为 suggested_params、reasoning、group_summary、sample_assessments。"
            " 重要：若在 reasoning 中提到任何参数调整，必须在 suggested_params 中给出对应 fastp 参数键与目标值，且仅包含需要变更的键（与基础参数相同的值不要包含）。"
            " sample_assessments 是一个字典，键为sample_id，值为质量评级（优秀/良好/合格/需改进/差）。"
        )

        # 构建提示词，避免f-string中的格式化问题
        base_params_str = _json.dumps(base_params, ensure_ascii=False)
        compact_str = _json.dumps(compact, ensure_ascii=False, indent=2)
        
        # 最近历史（最多3次）可选上下文；若未提供则为空
        history = history or []
        recent = history[-3:] if len(history) > 3 else history
        history_compact = []
        for history_entry in recent:
            history_compact.append({
                "version": history_entry.get("version"),
                "applied_params": history_entry.get("params", {}),
                "suggested": history_entry.get("optimized_params", {}),
                "success_rate": history_entry.get("execution_result", {}).get("success_rate", 0),
            })
        history_str = _json.dumps(history_compact, ensure_ascii=False, indent=2)
        
        param_prompt = f"""
基础参数（可作为默认或下限约束）：
{base_params_str}

样本统计（包含详细质量指标）：
{compact_str}

最近历史（请根据历史来决定参数，避免反复调节同样的参数）：
{history_str}

请基于以下指标评估每个样本的质量：
- reads_retention: reads保留率（>80%优秀，70-80%良好，60-70%合格，50-60%需改进，<50%差）
- q30_rate: Q30碱基比例（>90%优秀，80-90%良好，70-80%合格，60-70%需改进，<60%差）
- gc_content: GC含量变化（变化<5%正常，5-10%轻微异常，>10%异常）
- duplication_rate: 重复率（<10%优秀，10-20%良好，20-30%合格，30-40%需改进，>40%差）
- adapter_trimmed_percent: 接头污染（<5%优秀，5-10%良好，10-20%合格，20-30%需改进，>30%差）

请：
1) 评估每个样本的质量等级，在 sample_assessments 中返回 {{"sample_id": "质量评级"}} 的映射
2) 以中文总结整体质量、潜在问题、样本间差异，写入 group_summary
3) 给出谨慎且通用的 fastp 参数建议（仅列出需要变更的键及目标值），写入 suggested_params
4) **重要：必须在 reasoning 中提供参数决策理由**：
   - 如果有参数建议，请解释为什么建议这些变更及其预期效果
   - 如果无需调整参数，请说明"当前参数配置适合本批样本，质量指标符合预期"
   - reasoning字段不能为空，至少要有30个中文字符的解释

输出规范（严格遵守）：
- 仅输出以下 fastp 键（白名单）：
  qualified_quality_phred, unqualified_percent_limit, n_base_limit, length_required,
  adapter_trimming, quality_filtering, length_filtering, phred64, reads_to_process, fix_mgi_id,
  detect_adapter_for_pe, trim_front1, trim_tail1, max_len1, trim_front2, trim_tail2, max_len2,
  trim_poly_g, poly_g_min_len, disable_trim_poly_g, trim_poly_x, poly_x_min_len,
  cut_front, cut_tail, cut_right, cut_window_size, cut_mean_quality, cut_front_window_size,
  cut_front_mean_quality, cut_tail_window_size, cut_tail_mean_quality, cut_right_window_size,
  cut_right_mean_quality, average_qual, disable_length_filtering, length_limit, low_complexity_filter,
  complexity_threshold, correction, overlap_len_require, overlap_diff_limit, overlap_diff_percent_limit,
  overrepresentation_analysis, overrepresentation_sampling。
- 类型必须与基础参数一致：布尔值使用 true/false；整数/数值不得使用字符串表示；不得输出 null。
- suggested_params 只能包含“与基础参数 base_params 不同”的键值；不要回显相同值。
- 若 reasoning 出现“建议/启用/降低/提高/阈值/polyG/低复杂度”等字样，则 suggested_params 不得为空。
- 若确实无需调整，请返回空对象 {{}} 并在 reasoning 明确写“无需调整”。

示例（仅示意，不要照抄）：
{{
  "suggested_params": {{
    "trim_poly_g": true,
    "poly_g_min_len": 10,
    "low_complexity_filter": true,
    "complexity_threshold": 20
  }},
  "reasoning": "重复率偏高且末端质量下降，启用polyG并设置最小长度10；开启低复杂度过滤并降低阈值至20，尽量不影响reads保留率。",
  "group_summary": "总体质量优秀但存在重复率偏高…",
  "sample_assessments": {{"S1": "良好", "S2": "良好"}}
}}

最终仅输出 JSON 对象：{{"suggested_params":{{}},"reasoning":"","group_summary":"","sample_assessments":{{}}}}
"""
        
        try:
            llm = self.llm_manager.get_llm()
            resp = llm.invoke([HumanMessage(content=sys_prompt + "\n\n" + param_prompt)])
            content = resp.content if hasattr(resp, "content") else str(resp)
            
            # 解析JSON响应
            parsed = self._parse_llm_json_response(str(content))
            
            # 兜底结构校验
            if not isinstance(parsed, dict):
                raise ValueError("LLM输出非字典")
            parsed.setdefault("suggested_params", {})
            parsed.setdefault("reasoning", "")
            parsed.setdefault("group_summary", "")
            parsed.setdefault("sample_assessments", {})
            
            # 确保 reasoning 字段不为空
            if not parsed["reasoning"].strip():
                if parsed.get("suggested_params"):
                    parsed["reasoning"] = "基于质量指标分析，提供了参数优化建议以改善数据质量。"
                else:
                    parsed["reasoning"] = "当前参数配置适合本批样本，质量指标符合预期，无需调整。"
            
            return parsed
            
        except Exception as e:
            print(f"⚠️ LLM组级分析失败: {str(e)}")
            # 返回空结果，让调用方使用备选方法
            return {
                "suggested_params": {},
                "reasoning": f"LLM分析失败: {str(e)}",
                "group_summary": "",
                "sample_assessments": {}
            }
    
    
    def _get_base_params_from_nf(self, nextflow_config_values: Dict[str, Any]) -> Dict[str, Any]:
        """从 nextflow_config 提取（或默认生成）fastp 全量基础参数模板。

        说明：保留项目先前既定默认（如 Q20=20、length_required=15），
        其余参数补齐并给出与 fastp 语义一致的保守默认值。
        """
        return {
            # 核心质量/长度过滤（保持现有默认）
            "qualified_quality_phred": nextflow_config_values.get("qualified_quality_phred", 20),  # -q
            "unqualified_percent_limit": nextflow_config_values.get("unqualified_percent_limit", 40),  # -u
            "n_base_limit": nextflow_config_values.get("n_base_limit", 5),  # -n
            "length_required": nextflow_config_values.get("length_required", 15),  # -l
            "adapter_trimming": nextflow_config_values.get("adapter_trimming", True),
            "quality_filtering": nextflow_config_values.get("quality_filtering", True),
            "length_filtering": nextflow_config_values.get("length_filtering", True),  # 反向对应 -L

            # 输入质量与读取控制
            "phred64": nextflow_config_values.get("phred64", False),  # -6/--phred64
            "reads_to_process": nextflow_config_values.get("reads_to_process", 0),  # --reads_to_process (0=全部)
            "fix_mgi_id": nextflow_config_values.get("fix_mgi_id", False),  # --fix_mgi_id
            "detect_adapter_for_pe": nextflow_config_values.get("detect_adapter_for_pe", None),  # --detect_adapter_for_pe（默认None，批量逻辑另行决定）

            # 前后端定长修剪与最大长度
            "trim_front1": nextflow_config_values.get("trim_front1", 0),  # -f
            "trim_tail1": nextflow_config_values.get("trim_tail1", 0),  # -t
            "max_len1": nextflow_config_values.get("max_len1", 0),  # -b (0=无限制)
            "trim_front2": nextflow_config_values.get("trim_front2", 0),  # -F
            "trim_tail2": nextflow_config_values.get("trim_tail2", 0),  # -T
            "max_len2": nextflow_config_values.get("max_len2", 0),  # -B (0=无限制)

            # polyG / polyX 修剪
            "trim_poly_g": nextflow_config_values.get("trim_poly_g", False),  # -g
            "poly_g_min_len": nextflow_config_values.get("poly_g_min_len", None),  # --poly_g_min_len（仅在启用 -g 时生效）
            "disable_trim_poly_g": nextflow_config_values.get("disable_trim_poly_g", False),  # -G
            "trim_poly_x": nextflow_config_values.get("trim_poly_x", False),  # -x
            "poly_x_min_len": nextflow_config_values.get("poly_x_min_len", None),  # --poly_x_min_len

            # 滑窗切除开关（cut_* 开关使用布尔值控制启用/禁用）
            "cut_front": nextflow_config_values.get("cut_front", False),   # -5 启用前端切除
            "cut_tail": nextflow_config_values.get("cut_tail", False),     # -3 启用后端切除  
            "cut_right": nextflow_config_values.get("cut_right", False),   # -r 启用右端切除
            "cut_window_size": nextflow_config_values.get("cut_window_size", 4),  # -W
            "cut_mean_quality": nextflow_config_values.get("cut_mean_quality", 20),  # -M
            # 具体窗口与质量（未指定时沿用全局窗口与质量门限，给出显式默认以便模板完整）
            "cut_front_window_size": nextflow_config_values.get("cut_front_window_size", 4),
            "cut_front_mean_quality": nextflow_config_values.get("cut_front_mean_quality", 20),
            "cut_tail_window_size": nextflow_config_values.get("cut_tail_window_size", 4),
            "cut_tail_mean_quality": nextflow_config_values.get("cut_tail_mean_quality", 20),
            "cut_right_window_size": nextflow_config_values.get("cut_right_window_size", 4),
            "cut_right_mean_quality": nextflow_config_values.get("cut_right_mean_quality", 20),

            # 质量/长度过滤细化
            "average_qual": nextflow_config_values.get("average_qual", 0),  # -e (0=不启用)
            "disable_length_filtering": nextflow_config_values.get("disable_length_filtering", False),  # -L
            "length_limit": nextflow_config_values.get("length_limit", 0),  # --length_limit (0=无限制)
            "low_complexity_filter": nextflow_config_values.get("low_complexity_filter", False),  # -y
            "complexity_threshold": nextflow_config_values.get("complexity_threshold", 30),  # -Y

            # PE 重叠校正与检测
            "correction": nextflow_config_values.get("correction", False),  # -c（仅PE）
            "overlap_len_require": nextflow_config_values.get("overlap_len_require", 30),
            "overlap_diff_limit": nextflow_config_values.get("overlap_diff_limit", 5),
            "overlap_diff_percent_limit": nextflow_config_values.get("overlap_diff_percent_limit", 20),

            # 过表达序列分析
            "overrepresentation_analysis": nextflow_config_values.get("overrepresentation_analysis", False),  # -p
            "overrepresentation_sampling": nextflow_config_values.get("overrepresentation_sampling", 20),  # -P
        }


    def _summarize_batch(self, results: List[Dict[str, Any]], base_params: Dict[str, Any], optimized_params: Dict[str, Any], sample_assessments: Optional[Dict[str, str]] = None) -> str:
        """生成批次处理的总结报告"""
        total = len(results)
        successful = sum(1 for r in results if r.get("success"))
        summary = [
            "\n🧬 FastP质控处理完成",
            "",
            "📊 处理统计:",
            f"- 总样本数: {total}",
            f"- 成功处理: {successful}",
            f"- 失败样本: {total - successful}",
            "",
            "📋 样本详情:",
        ]

        for r in results:
            if r.get("success"):
                stats = r.get("statistics", {})
                retention = stats.get("reads_retention", 0) * 100
                q30_rate = stats.get("q30_rate", 0) * 100
                quality = sample_assessments.get(r['sample_id'], '未评级') if sample_assessments else '未评级'
                summary.append(f"- {r['sample_id']}: ✅ 成功 ({quality})")
                summary.append(f"  保留率: {retention:.1f}% | Q30: {q30_rate:.1f}%")
            else:
                summary.append(f"- {r['sample_id']}: ❌ 失败 - {r.get('error', 'unknown')}")

        summary.append("")
        summary.append("🧠 组级参数建议 (Fastp)：")
        changed = False
        for k in sorted(optimized_params.keys()):
            b = base_params.get(k)
            o = optimized_params.get(k)
            if b != o:
                summary.append(f"- {k}: {b} -> {o}")
                changed = True
        if not changed:
            summary.append("- 无需调整，保留默认参数")

        return "\n".join(summary)

    
    def run_batch(
        self,
        sample_groups: List[Dict[str, Any]],
        nextflow_config: Dict[str, Any],
        current_params: Optional[Dict[str, Any]] = None,
        version: int = 1,
        version_history: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """单次 Nextflow 批量执行 fastp，并生成组级建议与总结。
        
        Args:
            sample_groups: 样本组信息
            nextflow_config: Nextflow配置
            current_params: 当前运行参数（调用方传入，推荐来自 state.fastp_params）
            version: 参数版本号，用于文件命名
        """
        # 存储全局配置供其他方法使用
        self._global_nextflow_config = nextflow_config or {}

        from pathlib import Path

        # 结果根目录：使用 nextflow_config.results_dir
        results_root = Path(nextflow_config["results_dir"]) / "fastp"
        results_root.mkdir(parents=True, exist_ok=True)

        # 1. 优先使用current_params（来自state的历史优化参数）
        # 2. 其次使用nextflow_config中的显式配置
        # 3. 最后使用系统默认参数
        base_params = self._get_base_params_from_nf(nextflow_config or {})
        
        # 如果有历史优化参数，优先使用
        if current_params:
            print(f"🔄 使用历史优化参数: {len(current_params)} 个参数")
            effective_params = base_params.copy()
            effective_params.update(current_params)
        else:
            print(f"🆕 使用默认参数配置")
            effective_params = base_params

        # CPU资源：仅通过 Nextflow 的 fastp_cpus 管理
        fastp_cpus = (nextflow_config or {}).get('fastp_cpus')

        # detect_adapter_for_pe：优先使用 fastp_params（effective_params）；缺省时若任一组为PE则默认开启
        any_paired = any(bool(sample_group.get("read2")) for sample_group in (sample_groups or []))

        # 构建批量 params 字典
        batch_params: Dict[str, Any] = {
            "sample_groups": sample_groups or [],
            "results_dir": str(results_root.parent),  # fastp.nf 会拼接 fastp/${sample_id}
        }
        if fastp_cpus is not None:
            batch_params["fastp_cpus"] = int(fastp_cpus)

        # 合并有效 fastp 参数（优先使用历史优化参数）
        batch_params.update({
            "qualified_quality_phred": effective_params.get("qualified_quality_phred", 20),
            "unqualified_percent_limit": effective_params.get("unqualified_percent_limit", 40),
            "n_base_limit": effective_params.get("n_base_limit", 5),
            "length_required": effective_params.get("length_required", 15),
            "adapter_trimming": effective_params.get("adapter_trimming", True),
            "quality_filtering": effective_params.get("quality_filtering", True),
            "length_filtering": effective_params.get("length_filtering", True),
            "detect_adapter_for_pe": effective_params.get("detect_adapter_for_pe") if effective_params.get("detect_adapter_for_pe") is not None else any_paired,
        })

        # 添加高级参数（优先从effective_params读取）
        advanced_params = {
            "average_qual", "length_limit", "low_complexity_filter", "trim_poly_g", "disable_trim_poly_g",
            "trim_poly_x", "poly_g_min_len", "poly_x_min_len", "cut_front", "cut_tail",
            "cut_right", "cut_window_size", "cut_mean_quality", "cut_front_window_size",
            "cut_front_mean_quality", "cut_tail_window_size", "cut_tail_mean_quality",
            "cut_right_window_size", "cut_right_mean_quality", "correction", "complexity_threshold",
            "overlap_len_require", "overlap_diff_limit", "overlap_diff_percent_limit",
            "overrepresentation_sampling"
        }
        
        # 从effective_params中添加高级参数
        for param in advanced_params:
            if param in effective_params:
                batch_params[param] = effective_params[param]

        # 不再让 nextflow_config 覆盖 fastp 参数：fastp_params 是唯一真相源
        # 若需要通过 nextflow_config 显式注入 fastp 相关键，应在 Prepare 阶段同步写入 fastp_params

        # 写入版本化参数文件
        params_file = results_root / f"fastp_params.v{version}.json"
        params_file_latest = results_root / "fastp_params.batch.json"  # Nextflow执行用的文件
        
        print(f"📋 保存参数配置: v{version} (版本化管理)")
        
        # 保存版本化参数文件（带时间戳和版本信息）
        versioned_config = {
            "version": version,
            "timestamp": __import__('datetime').datetime.now().isoformat(),
            "description": f"FastP参数配置 v{version}",
            "config": batch_params.copy()
        }
        
        with open(params_file, 'w', encoding='utf-8') as f:
            json.dump(versioned_config, f, indent=2, ensure_ascii=False)
        
        # 同时保存为Nextflow执行用的标准文件名
        with open(params_file_latest, 'w', encoding='utf-8') as f:
            json.dump(batch_params, f, indent=2, ensure_ascii=False)

        # 运行 Nextflow（使用标准文件名）
        cmd = ["nextflow", "run", "/fastp.nf", "-params-file", str(params_file_latest), "-resume"]
        print("🚀 执行批量 fastp:", " ".join(cmd))
        completed_process = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

        if completed_process.returncode != 0:
            print(f"⚠️ Nextflow返回非零状态码: {completed_process.returncode}")
            print(f"   stderr: {completed_process.stderr[:200]}...")
            # 但仍继续检查输出文件，因为可能只是警告
        
        # 检查是否有样本成功产生输出文件
        successful_samples = 0
        for sample_group in sample_groups or []:
            sample_id = sample_group.get("sample_id", "unknown")
            sample_dir = results_root / sample_id
            json_file = sample_dir / f"{sample_id}.fastp.json"
            if json_file.exists():
                successful_samples += 1
        
        if successful_samples == 0:
            # 完全失败：返回统一错误
            return {
                "samples": [
                    {"sample_id": sample_group.get("sample_id", "unknown"), "success": False, "error": completed_process.stderr}
                    for sample_group in (sample_groups or [])
                ],
                "default_params": base_params,
                "optimized_params": {},
                "summary": f"批量执行失败: {completed_process.stderr[:200]}",
                "reasoning": "",
                "total": len(sample_groups or []),
                "success_count": 0,
                "version": version,
                "version_files": {
                    "versioned": str(params_file),
                    "latest": str(params_file_latest)
                },
            }
        elif completed_process.returncode != 0:
            print(f"🔧 Nextflow有警告但{successful_samples}/{len(sample_groups)}样本成功，继续分析")

        # 成功：逐样本收集报告并做分析
        results: List[Dict[str, Any]] = []
        for sample_group in sample_groups or []:
            sample_id = sample_group.get("sample_id") or "unknown"
            sample_output_dir = results_root / sample_id
            # 构造 FastpResult 并从报告填充
            fastp_result = FastpResult(success=True, exit_code=0)
            # 报告路径
            html_report_path = sample_output_dir / f"{sample_id}.fastp.html"
            json_report_path = sample_output_dir / f"{sample_id}.fastp.json"
            if html_report_path.exists():
                fastp_result.html_report = str(html_report_path)
            if json_report_path.exists():
                fastp_result.json_report = str(json_report_path)
                fastp_result = self._parse_json_report_with_path(str(json_report_path), fastp_result)
            else:
                fastp_result.success = False
                fastp_result.error_message = "缺少 fastp JSON 报告"

            # 输出文件
            if sample_group.get("read2"):
                # PE
                trimmed_read1_path = sample_output_dir / f"{sample_id}_1.trimmed.fastq.gz"
                trimmed_read2_path = sample_output_dir / f"{sample_id}_2.trimmed.fastq.gz"
                output_files_list: List[str] = []
                if trimmed_read1_path.exists():
                    output_files_list.append(str(trimmed_read1_path))
                if trimmed_read2_path.exists():
                    output_files_list.append(str(trimmed_read2_path))
                fastp_result.output_files = output_files_list
            else:
                # SE
                single_trimmed_path = sample_output_dir / f"{sample_id}.single.trimmed.fastq.gz"
                if single_trimmed_path.exists():
                    fastp_result.output_files = [str(single_trimmed_path)]


            sample_result: Dict[str, Any] = {
                "sample_id": sample_id,
                "success": fastp_result.success,
                "statistics": {
                    # 基本统计
                    "reads_before": fastp_result.total_reads_before,
                    "reads_after": fastp_result.total_reads_after,
                    "reads_retention": fastp_result.total_reads_after / fastp_result.total_reads_before if fastp_result.total_reads_before > 0 else 0,
                    "q30_rate": fastp_result.q30_bases_after / fastp_result.total_bases_after if fastp_result.total_bases_after > 0 else 0,
                    # 额外的质量指标
                    "gc_content_before": fastp_result.gc_content_before,
                    "gc_content_after": fastp_result.gc_content_after,
                    "mean_length_before": fastp_result.read_length_before,
                    "mean_length_after": fastp_result.read_length_after,
                    "duplication_rate": fastp_result.duplication_rate,
                    "adapter_trimmed_reads": fastp_result.adapter_trimmed_reads,
                    "adapter_trimmed_percent": fastp_result.adapter_trimmed_reads / fastp_result.total_reads_before if fastp_result.total_reads_before > 0 else 0,
                    # 过滤原因统计
                    "low_quality_reads": fastp_result.low_quality_reads,
                    "too_short_reads": fastp_result.too_short_reads,
                    "too_many_n_reads": fastp_result.too_many_n_reads,
                },
                "output_files": fastp_result.output_files,
                "reports": {"html": fastp_result.html_report, "json": fastp_result.json_report},
            }
            if not fastp_result.success:
                sample_result["error"] = fastp_result.error_message

            results.append(sample_result)

        # 优先：一次性 LLM 组级总结 + 参数建议
        try:
            # 简化：传入基础参数模板，不使用当前基线做对比
            group_llm = self._get_llm_group_summary_and_params(results, base_params, history=version_history)
            optimized_params = group_llm.get("suggested_params", {}) or {}
            reasoning = group_llm.get("reasoning", "") or ""
            summary = group_llm.get("group_summary") or ""
            sample_assessments = group_llm.get("sample_assessments", {}) or {}
            
            # 如果没有summary，使用默认总结
            if not summary:
                summary = self._summarize_batch(results, base_params, optimized_params, sample_assessments)
            
            # 如果LLM没有提供有效建议，使用默认参数
            if not optimized_params and not reasoning:
                optimized_params = {}
                reasoning = "LLM未提供参数建议，保持默认参数"
                
        except Exception as e:
            print(f"⚠️ 组级LLM总结/参数建议失败，使用回退逻辑: {str(e)}")
            # 回退：使用默认参数
            optimized_params = {}
            reasoning = "LLM分析失败，保持默认参数"
            summary = self._summarize_batch(results, effective_params, optimized_params, {})

        # 计算下次迭代的参数：直接合并 LLM 建议（不做实时差异过滤）
        next_params = effective_params.copy()
        if optimized_params:
            next_params.update(optimized_params)
            print(f"📈 参数进化: {len(optimized_params)} 个参数将更新")
        else:
            print(f"📊 参数稳定: 保持当前配置")
            # 无建议时保持原理由

        return {
            "samples": results,
            "default_params": base_params,            # 系统默认参数（不变）
            "current_params": effective_params,       # 本次执行使用的参数
            "optimized_params": optimized_params,     # 直接返回LLM建议（不做差异过滤）
            "next_params": next_params,               # 下次执行应使用的参数
            "version": version,  # 当前参数版本号
            "version_files": {  # 版本化文件路径
                "versioned": str(params_file),
                "latest": str(params_file_latest)
            },
            "summary": summary,
            "reasoning": reasoning,
            "total": len(results),
            "success_count": sum(1 for r in results if r.get('success')),
        }
