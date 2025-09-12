"""User Confirm 节点重构模块 - 纯数字选择模式

提供用户确认界面的结构化、可测试的实现，采用纯数字选择交互模式。

模块结构：
- view_model: 数据模型定义 (包含CommandHint的index字段)
- presenter: 视图模型构建逻辑 (支持数字索引分配)
- render: 渲染器，将模型转为字符串 (数字化渲染)
- commands: 命令解析器 (纯数字选择解析)
- diffing: 参数差异对比工具
"""

from .view_model import ConfirmView, Section, CommandHint, ConfirmDecision
from .presenter import build_confirm_view
from .render import render_confirm
from .commands import (
    parse_numeric_selection, parse_execution_mode_selection,
    get_execution_mode_selection, parse_confirm_command
)
from .diffing import flatten_dict, diff_base_mods_opt, ParamDiff

__all__ = [
    'ConfirmView', 'Section', 'ParamDiff', 'CommandHint', 'ConfirmDecision',
    'build_confirm_view',
    'render_confirm',
    'parse_numeric_selection', 'parse_execution_mode_selection', 
    'get_execution_mode_selection', 'parse_confirm_command',
    'flatten_dict', 'diff_base_mods_opt'
]