"""User Confirm 节点重构模块

提供用户确认界面的结构化、可测试的实现，降低复杂度并提升可维护性。

模块结构：
- view_model: 数据模型定义
- presenter: 视图模型构建逻辑
- render: 渲染器，将模型转为字符串
- commands: 命令解析器
- diffing: 参数差异对比工具
"""

from .view_model import ConfirmView, Section, CommandHint
from .presenter import build_confirm_view
from .render import render_confirm
from .commands import parse_confirm_command, ConfirmDecision, get_available_commands, validate_command
from .diffing import flatten_dict, diff_base_mods_opt, ParamDiff

__all__ = [
    'ConfirmView', 'Section', 'ParamDiff', 'CommandHint',
    'build_confirm_view',
    'render_confirm',
    'parse_confirm_command', 'ConfirmDecision', 'get_available_commands', 'validate_command',
    'flatten_dict', 'diff_base_mods_opt'
]