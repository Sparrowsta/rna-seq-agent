"""命令解析器 - 纯数字选择模式

处理用户确认界面的纯数字选择输入，完全弃用斜杠命令。
将数字选择转换为标准化的决策对象。
"""

from typing import Dict, Any, List, Optional
from .view_model import ConfirmDecision, CommandHint


def parse_numeric_selection(
    raw_input: str, 
    commands: List[CommandHint], 
    context: Dict[str, Any]
) -> ConfirmDecision:
    """
    解析用户的纯数字选择输入
    
    Args:
        raw_input: 原始用户输入
        commands: 可用命令列表 (带index)
        context: 上下文信息
        
    Returns:
        ConfirmDecision: 标准化的决策对象
    """
    # 规范化输入
    user_input = raw_input.strip()
    
    # 检查是否为纯数字
    try:
        selected_number = int(user_input)
    except ValueError:
        # 非数字输入，返回无效
        return ConfirmDecision(
            decision='cancel',
            payload={
                'error': 'invalid_input',
                'message': f'无效输入："{raw_input}"。请输入数字选择操作。',
                'raw_input': raw_input
            }
        )
    
    # 0表示返回/取消
    if selected_number == 0:
        return ConfirmDecision(
            decision='cancel',
            payload={'raw_input': raw_input, 'user_cancel': True}
        )
    
    # 查找对应的命令
    selected_command = None
    for cmd in commands:
        if cmd.available and cmd.index == selected_number:
            selected_command = cmd
            break
    
    if not selected_command:
        # 越界或无效数字
        available_numbers = [str(cmd.index) for cmd in commands if cmd.available and cmd.index]
        return ConfirmDecision(
            decision='cancel',
            payload={
                'error': 'out_of_range',
                'message': f'无效选择：{selected_number}。可用选项：{", ".join(available_numbers)}, 0(取消)',
                'raw_input': raw_input,
                'available_numbers': available_numbers
            }
        )
    
    # 将选中的命令转换为决策
    return _convert_command_to_decision(selected_command, context, raw_input)


def _convert_command_to_decision(
    command: CommandHint, 
    context: Dict[str, Any], 
    raw_input: str
) -> ConfirmDecision:
    """将CommandHint转换为ConfirmDecision"""
    
    completed_steps = context.get('completed_steps', [])
    current_step = context.get('current_step', '')
    
    # 根据命令类型转换
    if command.command == "/execute_opt":
        return ConfirmDecision(
            decision='execute',
            execution_mode=None,  # 需要进一步选择模式
            payload={
                'raw_input': raw_input,
                'needs_mode_selection': True,
                'command_description': command.description
            }
        )
    
    elif command.command == "/yolo":
        return ConfirmDecision(
            decision='execute',
            execution_mode='yolo',
            payload={
                'mode_description': '🎯 YOLO自动模式：全程自动执行分析流程',
                'raw_input': raw_input
            }
        )
    
    elif command.command == "/continue":
        # 根据完成步骤确定下一步
        if "featurecounts" in completed_steps:
            decision = 'continue_analysis'
        elif "star" in completed_steps:
            decision = 'continue_featurecounts'
        elif "fastp" in completed_steps:
            decision = 'continue_star'
        else:
            decision = 'continue'
        
        return ConfirmDecision(
            decision=decision,
            payload={'raw_input': raw_input}
        )
    
    elif command.command == "/restart":
        return ConfirmDecision(
            decision='execute',
            execution_mode='single',
            payload={'restart': True, 'raw_input': raw_input}
        )
    
    elif command.command == "/re_opt":
        # 确定目标步骤
        target_step = current_step if current_step in {"fastp", "star", "featurecounts"} else None
        if not target_step and completed_steps:
            # 从已完成步骤中选择最后一个
            for step in reversed(["featurecounts", "star", "fastp"]):
                if step in completed_steps:
                    target_step = step
                    break
        
        if not target_step:
            target_step = 'fastp'  # 默认
        
        return ConfirmDecision(
            decision=target_step,
            execution_mode='optimized',
            payload={
                're_optimization': True,
                'target_step': target_step,
                'raw_input': raw_input
            }
        )
    
    elif command.command == "/modify":
        return ConfirmDecision(
            decision='modify',
            payload={
                'raw_input': raw_input,
                'needs_modify_content': True  # 标记需要二次输入
            }
        )
    
    elif command.command == "/cancel":
        return ConfirmDecision(
            decision='cancel',
            payload={'raw_input': raw_input}
        )
    
    elif command.command == "/quit":
        return ConfirmDecision(
            decision='quit',
            payload={'raw_input': raw_input}
        )
    
    else:
        # 未知命令
        return ConfirmDecision(
            decision='cancel',
            payload={
                'error': 'unknown_command',
                'message': f'未知命令：{command.command}',
                'raw_input': raw_input
            }
        )


def get_execution_mode_selection() -> List[str]:
    """获取执行模式选择菜单"""
    return [
        "请选择执行模式 (1-3, 0返回):",
        "    1) Single - 单次执行",
        "    2) Optimized - 优化模式", 
        "    3) Batch Optimize - 批次优化",
        "    0) 返回上级菜单"
    ]


def parse_execution_mode_selection(user_input: str) -> Optional[str]:
    """解析执行模式选择"""
    user_input = user_input.strip()
    
    try:
        selected = int(user_input)
        if selected == 0:
            return 'cancel'
        elif selected == 1:
            return 'single'
        elif selected == 2:
            return 'optimized'
        elif selected == 3:
            return 'batch_optimize'
        else:
            return None
    except ValueError:
        return None


# 向后兼容函数 (如果其他地方还在调用)
def parse_confirm_command(raw_input: str, context: Dict[str, Any]) -> ConfirmDecision:
    """向后兼容的解析函数 - 重定向到错误处理"""
    # 使用context参数避免未使用警告
    _ = context
    return ConfirmDecision(
        decision='cancel',
        payload={
            'error': 'deprecated_function',
            'message': '旧版命令解析已弃用，请使用纯数字选择模式',
            'raw_input': raw_input
        }
    )