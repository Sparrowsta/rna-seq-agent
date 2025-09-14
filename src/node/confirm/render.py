"""渲染器

将结构化的视图模型转换为字符串列表，保持现有的中文风格和图标。
负责格式化输出，不包含任何业务逻辑。
"""

from typing import List
from .view_model import ConfirmView, Section, SummaryItem, ResourceItem, CommandHint, ParamItem


def render_confirm(view: ConfirmView) -> List[str]:
    """
    将ConfirmView渲染为字符串列表
    
    Args:
        view: 视图模型对象
        
    Returns:
        格式化的字符串列表，每个元素为一行输出
    """
    lines = []
    
    # 页面标题（轻量版样式，无首行空白）
    lines.extend(["-" * 60, "🎯 分析配置确认", "-" * 60])
    
    # 配置摘要
    lines.extend(_render_summary(view.summary))
    
    # 资源配置
    lines.extend(_render_resources(view.resources))
    
    # 配置对比（Nextflow和Resource的差异对比已在presenter中处理）
    # 这里主要渲染参数三层对比区域
    for section in view.sections:
        lines.extend(_render_section(section))
    
    # 配置理由
    if view.config_reasoning:
        lines.append("")
        lines.append("💭 配置理由:")
        lines.append(f"   {view.config_reasoning}")
    
    # 执行进度信息
    lines.extend(_render_progress(view))
    
    # 命令提示
    lines.extend(_render_commands(view.commands))
    
    # 结尾细分割线
    lines.append("-" * 60)
    
    return lines


def _render_summary(summary_items: List[SummaryItem]) -> List[str]:
    """渲染配置摘要"""
    if not summary_items:
        return ["", "📋 配置摘要: ⚠️ 无配置信息"]
    
    lines = ["", "📋 配置摘要"]
    
    for item in summary_items:
        if not item.visible:
            continue
            
        # 特殊处理样本文件的显示
        if item.key == 'sample_groups' and isinstance(item.value, str) and item.value.endswith('个样本'):
            lines.append(f"   {item.icon} {item.label}: {item.value}")
            # 这里不能展示具体样本详情，因为没有原始数据
        else:
            lines.append(f"   {item.icon} {item.label}: {item.value}")
    
    return lines


def _render_resources(resource_items: List[ResourceItem]) -> List[str]:
    """渲染资源配置"""
    if not resource_items:
        return ["", "🖥️ 资源配置: 使用默认设置"]
    
    lines = ["", "🖥️ 资源配置"]
    
    for item in resource_items:
        lines.append(f"   {item.display_name}: {item.cpus}核, {item.memory}")
        if item.reasoning:
            lines.append(f"      💭 {item.reasoning}")
    
    return lines


def _render_section(section: Section) -> List[str]:
    """渲染参数对比区域"""
    if not section.visible:
        return []
    
    lines = ["", f"{section.icon} {section.title}"]
    
    # Effective（当前生效）
    if section.effective:
        lines.append(f"   📋 当前（Effective）:")
        for item in section.effective:
            lines.append(f"     - {item.key}: {item.value}")
    
    # User Mods（用户修改）
    if section.user_mods:
        lines.append("")
        lines.append(f"   ✏️ 用户修改（Mods）：")
        for item in section.user_mods:
            # 这里需要显示变化，但我们在ParamItem中没有old_value
            # 暂时只显示当前值
            lines.append(f"     - {item.key}: {item.value}")
    else:
        lines.append("")
        lines.append(f"   ✏️ 用户修改（Mods）： 无")
    
    # Optimizations（优化建议）
    if section.optimizations:
        lines.append("")
        lines.append(f"   ⚙️ 优化建议（Opt）：")
        for item in section.optimizations:
            if getattr(item, 'applied_optimization', False):
                suffix = " (已应用优化)"
            else:
                suffix = ""
            lines.append(f"     - {item.key}: {item.value}{suffix}")
    else:
        lines.append("")
        lines.append(f"   ⚙️ 优化建议（Opt）： 无")
    
    # 优化理由
    if section.reasoning_text:
        lines.append("")
        lines.append(f"   📝 优化理由：")
        for line in section.reasoning_text.splitlines():
            if line.strip():
                lines.append(f"     {line.strip()}")
    
    return lines


def _render_progress(view: ConfirmView) -> List[str]:
    """渲染执行进度信息"""
    lines = []
    
    # 批次优化完成提示
    if view.batch_optimization_complete and view.batch_optimizations_count > 0:
        lines.append("")
        lines.append(f"✅ 批次优化完成: 已收集{view.batch_optimizations_count}个工具的优化参数并应用")
    
    # 执行进度
    if view.completed_steps:
        lines.append("")
        lines.append(f"📊 执行进度: {' -> '.join(view.completed_steps)}")
        
        if view.current_step:
            lines.append(f"   🔄 当前步骤: {view.current_step}")
    
    return lines


def _render_commands(commands: List[CommandHint]) -> List[str]:
    """渲染命令提示"""
    if not commands:
        return []
    
    lines = ["", "🔄 请选择下一步操作:"]
    
    for cmd in commands:
        if cmd.available and cmd.index is not None:
            # 使用数字索引格式渲染
            icon_part = f" {cmd.icon}" if cmd.icon else ""
            lines.append(f"    {cmd.index}) {cmd.description}{icon_part}")
    
    return lines
