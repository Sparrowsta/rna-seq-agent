"""视图模型构建器

负责从 AgentState 构建结构化的 ConfirmView 对象，
聚合配置信息、参数对比、命令提示等，不包含任何 I/O 操作。
"""

from typing import Dict, Any, List, Optional
from ...state import AgentState
from .view_model import (
    ConfirmView, Section, SummaryItem, ResourceItem, 
    CommandHint, ParamItem
)
from .diffing import (
    diff_base_mods_opt, extract_user_modifications_from_history,
    extract_optimization_from_history, filter_ignored_keys
)


# 配置摘要映射表：键 -> (标签, 图标, 可见性谓词)
SUMMARY_MAPPING = {
    'genome_version': ('基因组版本', '🧬', lambda cfg: True),
    'species': ('物种', '🔬', lambda cfg: True),
    'qc_tool': ('质控工具', '🧹', lambda cfg: True),
    'align_tool': ('比对工具', '🎯', lambda cfg: True),
    'quant_tool': ('定量工具', '📊', lambda cfg: True),
    'paired_end': ('测序类型', '🔄', lambda cfg: True),
    'run_download_genome': ('下载基因组', '⬇️', lambda cfg: True),
    'run_build_star_index': ('构建STAR索引', '🏗️', lambda cfg: cfg.get('align_tool', '').lower() == 'star'),
    'run_build_hisat2_index': ('构建HISAT2索引', '🏗️', lambda cfg: cfg.get('align_tool', '').lower() == 'hisat2'),
    'sample_groups': ('样本文件', '📂', lambda cfg: True),
}

# 资源配置进程名称映射
RESOURCE_DISPLAY_MAPPING = {
    'prepare_star_index': ('STAR索引构建', '🏗️'),
    'prepare_hisat2_index': ('HISAT2索引构建', '🏗️'),
    'run_alignment': ('序列比对', '🎯'),
    'run_quality_control': ('质控处理', '🧹'),
    'run_quantification': ('基因定量', '📊'),
    'download_genome_fasta': ('FASTA下载', '⬇️'),
    'download_genome_gtf': ('GTF下载', '⬇️'),
}


def build_confirm_view(state: AgentState) -> ConfirmView:
    """
    从 AgentState 构建完整的用户确认视图模型
    
    Args:
        state: Agent状态对象
        
    Returns:
        ConfirmView: 结构化的视图模型
    """
    # 获取基础配置信息
    nextflow_config = state.nextflow_config or {}
    resource_config = state.resource_config or {}
    config_reasoning = getattr(state, 'config_reasoning', '') or ''
    
    # 构建各个部分
    summary = _build_summary(nextflow_config)
    resources = _build_resources(resource_config)
    sections = _build_param_sections(state, nextflow_config)
    commands = _build_commands(state)
    
    # 执行进度信息
    completed_steps = getattr(state, 'completed_steps', []) or []
    current_step = getattr(state, 'current_step', None)
    
    # 批次优化状态
    batch_optimizations = getattr(state, 'batch_optimizations', {}) or {}
    batch_optimization_complete = getattr(state, 'batch_optimization_complete', False)
    
    return ConfirmView(
        summary=summary,
        resources=resources,
        sections=sections,
        config_reasoning=config_reasoning.strip() or None,
        commands=commands,
        completed_steps=completed_steps,
        current_step=current_step,
        batch_optimization_complete=batch_optimization_complete,
        batch_optimizations_count=len(batch_optimizations)
    )


def _build_summary(nextflow_config: Dict[str, Any]) -> List[SummaryItem]:
    """构建配置摘要列表"""
    summary_items = []
    
    for key, (label, icon, visible_if) in SUMMARY_MAPPING.items():
        if not visible_if(nextflow_config):
            continue
            
        value = nextflow_config.get(key)
        if value is None:
            continue
            
        # 特殊处理某些值的显示格式
        display_value = _format_summary_value(key, value)
        
        summary_items.append(SummaryItem(
            key=key,
            value=display_value,
            label=label,
            icon=icon,
            visible=True
        ))
    
    return summary_items


def _format_summary_value(key: str, value: Any) -> Any:
    """格式化摘要值的显示"""
    if key == 'paired_end':
        return "双端测序" if value else "单端测序"
    elif key in ['run_download_genome', 'run_build_star_index', 'run_build_hisat2_index']:
        return "是" if value else "否"
    elif key == 'sample_groups' and isinstance(value, list):
        return f"{len(value)}个样本"
    return value


def _build_resources(resource_config: Dict[str, Dict[str, Any]]) -> List[ResourceItem]:
    """构建资源配置列表"""
    resource_items = []
    
    for process_name, config in resource_config.items():
        display_name, icon = RESOURCE_DISPLAY_MAPPING.get(process_name, (process_name, '⚙️'))
        
        cpus = config.get('cpus', 'N/A')
        memory = config.get('memory', 'N/A')
        reasoning = config.get('reasoning', '')
        
        resource_items.append(ResourceItem(
            process_name=process_name,
            display_name=f"{icon} {display_name}",
            cpus=str(cpus),
            memory=str(memory),
            reasoning=reasoning.strip() or None
        ))
    
    return resource_items


def _build_param_sections(state: AgentState, nextflow_config: Dict[str, Any]) -> List[Section]:
    """构建参数对比区域列表"""
    sections = []
    
    # FastP 参数对比（容错：若未声明qc_tool但存在参数同样显示）
    qc_tool = (nextflow_config.get('qc_tool') or '').lower()
    if qc_tool == 'fastp' or (getattr(state, 'fastp_params', {}) or {}):
        fastp_section = _build_fastp_section(state)
        if fastp_section:
            sections.append(fastp_section)
    
    # 比对工具参数对比：严格仅展示当前选择的比对器，杜绝同时显示
    align_tool = (nextflow_config.get('align_tool') or '').lower()

    if align_tool == 'hisat2':
        hisat2_section = _build_hisat2_section(state)
        if hisat2_section:
            sections.append(hisat2_section)
    elif align_tool == 'star':
        star_section = _build_star_section(state)
        if star_section:
            sections.append(star_section)
    # 未指定 align_tool 时，不展示任何比对器参数区，避免误导
    
    # FeatureCounts 参数对比（容错）
    quant_tool = (nextflow_config.get('quant_tool') or '').lower()
    if quant_tool == 'featurecounts' or (getattr(state, 'featurecounts_params', {}) or {}):
        fc_section = _build_featurecounts_section(state)
        if fc_section:
            sections.append(fc_section)
    
    return sections


def _build_fastp_section(state: AgentState) -> Optional[Section]:
    """构建FastP参数对比区域"""
    effective_params = getattr(state, 'fastp_params', {}) or {}
    if not effective_params:
        return None
    
    base_params = getattr(state, 'prepare_defaults_fastp_params', {}) or {}
    if not base_params:
        base_params = dict(effective_params)
    
    # 用户修改记录
    modification_history = getattr(state, 'modification_history', []) or []
    user_mods = extract_user_modifications_from_history(modification_history, 'fastp')
    
    # 优化记录：优先使用新字段 fastp_optimization_params；否则回退到历史
    opt_applied = getattr(state, 'fastp_optimization_params', None)
    if not opt_applied:
        fastp_history = getattr(state, 'fastp_params_history', []) or []
        opt_applied = extract_optimization_from_history(fastp_history)
    
    # 优化理由文本
    reasoning_text = getattr(state, 'fastp_optimization_suggestions', '') or ''
    
    # 构建三层对比
    param_diff = diff_base_mods_opt(base_params, effective_params, user_mods, opt_applied)
    
    return Section(
        title="FastP 参数（三层对比）",
        icon="🧹",
        effective=_create_param_items(param_diff.effective, "当前生效"),
        user_mods=_create_param_items(param_diff.user_modifications, "用户修改") if param_diff.user_modifications else [],
        optimizations=_create_param_items(param_diff.optimizations, "优化建议") if param_diff.optimizations else [],
        reasoning_text=reasoning_text.strip() or None,
        visible=True
    )


def _build_star_section(state: AgentState) -> Optional[Section]:
    """构建STAR参数对比区域"""
    effective_params = getattr(state, 'star_params', {}) or {}
    if not effective_params:
        return None
    
    base_params = getattr(state, 'prepare_defaults_star_params', {}) or {}
    if not base_params:
        base_params = dict(effective_params)
    
    # 用户修改记录
    modification_history = getattr(state, 'modification_history', []) or []
    user_mods = extract_user_modifications_from_history(modification_history, 'star')
    
    # 优化记录：优先使用新字段 star_optimization_params；否则回退到历史
    opt_applied = getattr(state, 'star_optimization_params', None)
    if not opt_applied:
        star_history = getattr(state, 'star_params_history', []) or []
        opt_applied = extract_optimization_from_history(star_history)
    
    # 优化理由文本
    reasoning_text = getattr(state, 'star_optimization_suggestions', '') or ''
    
    # 构建三层对比
    param_diff = diff_base_mods_opt(base_params, effective_params, user_mods, opt_applied)
    
    return Section(
        title="STAR 参数（三层对比）",
        icon="🎯",
        effective=_create_param_items(param_diff.effective, "当前生效"),
        user_mods=_create_param_items(param_diff.user_modifications, "用户修改") if param_diff.user_modifications else [],
        optimizations=_create_param_items(param_diff.optimizations, "优化建议") if param_diff.optimizations else [],
        reasoning_text=reasoning_text.strip() or None,
        visible=True
    )


def _build_hisat2_section(state: AgentState) -> Optional[Section]:
    """构建HISAT2参数对比区域"""
    effective_params = getattr(state, 'hisat2_params', {}) or {}
    if not effective_params:
        return None

    base_params = getattr(state, 'prepare_defaults_hisat2_params', {}) or {}
    if not base_params:
        base_params = dict(effective_params)

    # 用户修改记录
    modification_history = getattr(state, 'modification_history', []) or []
    user_mods = extract_user_modifications_from_history(modification_history, 'hisat2')

    # 优化记录：优先使用新字段 hisat2_optimization_params；否则回退到历史（若存在）
    opt_applied = getattr(state, 'hisat2_optimization_params', None)
    if not opt_applied:
        hisat2_history = getattr(state, 'hisat2_params_history', []) or []
        opt_applied = extract_optimization_from_history(hisat2_history)

    # 优化理由文本
    reasoning_text = getattr(state, 'hisat2_optimization_suggestions', '') or ''

    # 构建三层对比
    param_diff = diff_base_mods_opt(base_params, effective_params, user_mods, opt_applied)

    return Section(
        title="HISAT2 参数（三层对比）",
        icon="🎯",
        effective=_create_param_items(param_diff.effective, "当前生效"),
        user_mods=_create_param_items(param_diff.user_modifications, "用户修改") if param_diff.user_modifications else [],
        optimizations=_create_param_items(param_diff.optimizations, "优化建议") if param_diff.optimizations else [],
        reasoning_text=reasoning_text.strip() or None,
        visible=True
    )


def _build_featurecounts_section(state: AgentState) -> Optional[Section]:
    """构建FeatureCounts参数对比区域"""
    effective_params = getattr(state, 'featurecounts_params', {}) or {}
    if not effective_params:
        return None
    
    base_params = getattr(state, 'prepare_defaults_featurecounts_params', {}) or {}
    if not base_params:
        base_params = dict(effective_params)
    
    # 用户修改记录
    modification_history = getattr(state, 'modification_history', []) or []
    user_mods = extract_user_modifications_from_history(modification_history, 'featurecounts')
    
    # 优化记录：优先使用新字段 featurecounts_optimization_params；否则回退到历史
    opt_applied = getattr(state, 'featurecounts_optimization_params', None)
    if not opt_applied:
        fc_history = getattr(state, 'featurecounts_params_history', []) or []
        opt_applied = extract_optimization_from_history(fc_history)
    
    # 优化理由文本
    reasoning_text = getattr(state, 'featurecounts_optimization_suggestions', '') or ''
    
    # 构建三层对比
    param_diff = diff_base_mods_opt(base_params, effective_params, user_mods, opt_applied)
    
    return Section(
        title="FeatureCounts 参数（三层对比）",
        icon="📊",
        effective=_create_param_items(param_diff.effective, "当前生效"),
        user_mods=_create_param_items(param_diff.user_modifications, "用户修改") if param_diff.user_modifications else [],
        optimizations=_create_param_items(param_diff.optimizations, "优化建议") if param_diff.optimizations else [],
        reasoning_text=reasoning_text.strip() or None,
        visible=True
    )


def _create_param_items(params: Dict[str, Any], category: str) -> List[ParamItem]:
    """从参数字典创建ParamItem列表"""
    items = []
    for key in sorted(params.keys()):
        items.append(ParamItem(
            key=key,
            value=params[key],
            applied_optimization=category == "优化建议"
        ))
    return items


def _build_commands(state: AgentState) -> List[CommandHint]:
    """构建动态命令提示列表"""
    commands = []
    
    completed_steps = getattr(state, 'completed_steps', []) or []
    current_step = getattr(state, 'current_step', '')
    
    # 根据执行进度动态生成命令
    if completed_steps:
        # 有执行进度，显示continue选项
        if "featurecounts" in completed_steps:
            commands.append(CommandHint(
                command="/continue", 
                description="继续到综合分析",
                icon="➡️"
            ))
        elif "star" in completed_steps:
            commands.append(CommandHint(
                command="/continue",
                description="继续到FeatureCounts定量",
                icon="➡️"
            ))
        elif "fastp" in completed_steps:
            commands.append(CommandHint(
                command="/continue",
                description="继续到STAR比对", 
                icon="➡️"
            ))
        
        commands.append(CommandHint(
            command="/restart",
            description="重新开始完整流水线",
            icon="🔄"
        ))
    else:
        # 无执行进度，显示execute选项
        commands.append(CommandHint(
            command="/execute_opt",
            description="执行分析（支持选择执行模式）",
            icon="⚡"
        ))
        # 添加YOLO自动模式选项
        commands.append(CommandHint(
            command="/yolo",
            description="YOLO 自动执行",
            icon="🎯"
        ))
    
    # 二次优化选项
    if current_step in {"fastp", "star", "featurecounts"}:
        commands.append(CommandHint(
            command="/re_opt",
            description=f"重新优化当前步骤({current_step})",
            icon="♻️"
        ))
    
    # 通用命令
    commands.extend([
        CommandHint(command="/modify", description="修改配置", icon="🔧"),
        CommandHint(command="/cancel", description="取消分析", icon="❌"),
        CommandHint(command="/quit", description="退出程序", icon="🚪"),
    ])
    
    # 为所有可用命令分配数字索引 (从1开始)
    index = 1
    for cmd in commands:
        if cmd.available:
            cmd.index = index
            index += 1
    
    return commands
