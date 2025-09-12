"""命令解析器

统一处理用户确认界面的命令输入，支持等价命令和智能解析。
将用户输入转换为标准化的决策对象。
"""

from typing import Dict, Any, List, Optional
from .view_model import ConfirmDecision


# 命令等价表：命令变体 -> 标准命令
COMMAND_ALIASES = {
    # execute 系列
    '/execute_opt': 'execute',
    '/execute': 'execute', 
    '/执行': 'execute',
    '/运行': 'execute',
    
    # YOLO 系列 - 新增
    '/yolo': 'yolo',
    '/自动': 'yolo',
    '/全自动': 'yolo',
    
    # continue 系列
    '/continue': 'continue',
    '/继续': 'continue',
    
    # restart 系列
    '/restart': 'restart',
    '/重启': 'restart',
    '/重新开始': 'restart',
    
    # re_opt 系列
    '/re_opt': 're_opt',
    '/重新优化': 're_opt',
    '/二次优化': 're_opt',
    
    # modify 系列  
    '/modify': 'modify',
    '/修改': 'modify',
    '/调整': 'modify',
    
    # cancel 系列
    '/cancel': 'cancel',
    '/取消': 'cancel',
    
    # quit 系列
    '/quit': 'quit',
    '/exit': 'quit',
    '/退出': 'quit',
    '/bye': 'quit',
}

# modify 系列前缀（用于检测带参数的modify命令）
MODIFY_PREFIXES = ['/modify', '/修改', '/调整']


def parse_confirm_command(raw_input: str, context: Dict[str, Any]) -> ConfirmDecision:
    """
    解析用户确认界面的命令输入
    
    Args:
        raw_input: 原始用户输入
        context: 上下文信息（包含completed_steps, current_step等）
        
    Returns:
        ConfirmDecision: 标准化的决策对象
    """
    # 规范化输入
    user_input = raw_input.strip()
    user_input_lower = user_input.lower()
    
    # 检测modify命令（可能带参数）
    modify_result = _parse_modify_command(user_input, user_input_lower)
    if modify_result:
        return modify_result
    
    # 检测标准命令
    standard_cmd = COMMAND_ALIASES.get(user_input_lower)
    if not standard_cmd:
        # 无效命令
        return ConfirmDecision(
            decision='cancel',  # 默认返回cancel，由调用方处理
            payload={'error': 'invalid_command', 'raw_input': raw_input}
        )
    
    # 根据标准命令和上下文构建决策
    return _build_decision(standard_cmd, context, raw_input)


def _parse_modify_command(user_input: str, user_input_lower: str) -> Optional[ConfirmDecision]:
    """解析modify命令，提取修改内容"""
    # 检测是否为modify命令
    is_modify = False
    modify_content = ""
    
    # 检查精确匹配
    if user_input_lower in MODIFY_PREFIXES:
        is_modify = True
        modify_content = ""
    else:
        # 检查前缀匹配（带参数）
        for prefix in MODIFY_PREFIXES:
            if user_input_lower.startswith(f"{prefix} "):
                is_modify = True
                modify_content = user_input[len(prefix):].strip()
                break
    
    if not is_modify:
        return None
    
    return ConfirmDecision(
        decision='modify',
        modify_content=modify_content or None,
        payload={
            'raw_input': modify_content or "",
            'has_content': bool(modify_content)
        }
    )


def _build_decision(
    command: str, 
    context: Dict[str, Any], 
    raw_input: str
) -> ConfirmDecision:
    """根据标准命令和上下文构建决策对象"""
    
    completed_steps = context.get('completed_steps', [])
    current_step = context.get('current_step', '')
    
    if command == 'execute':
        return _handle_execute_command(context, raw_input)
    
    elif command == 'yolo':
        # YOLO自动模式 - 新增
        return ConfirmDecision(
            decision='execute',
            execution_mode='yolo',
            payload={'mode_description': '🎯 YOLO自动模式：全程自动执行分析流程', 'raw_input': raw_input}
        )
    
    elif command == 'continue':
        return _handle_continue_command(completed_steps, current_step)
    
    elif command == 'restart':
        return ConfirmDecision(
            decision='execute',
            execution_mode='single',
            payload={'restart': True, 'raw_input': raw_input}
        )
    
    elif command == 're_opt':
        return _handle_re_opt_command(completed_steps, current_step, raw_input)
    
    elif command in ['cancel', 'quit']:
        return ConfirmDecision(
            decision=command,
            payload={'raw_input': raw_input}
        )
    
    else:
        # 未知命令，返回cancel
        return ConfirmDecision(
            decision='cancel',
            payload={'error': 'unknown_command', 'raw_input': raw_input}
        )


def _handle_execute_command(context: Dict[str, Any], raw_input: str) -> ConfirmDecision:
    """
    处理execute命令，可能需要弹出模式选择
    这里返回基础决策，具体的模式选择由调用方处理
    """
    # 使用context参数避免未使用警告
    _ = context  
    return ConfirmDecision(
        decision='execute',
        execution_mode=None,  # 需要进一步选择模式
        payload={'raw_input': raw_input, 'needs_mode_selection': True}
    )


def _handle_continue_command(completed_steps: List[str], current_step: str) -> ConfirmDecision:
    """处理continue命令，根据进度确定具体的continue类型"""
    # 使用current_step参数避免未使用警告
    _ = current_step
    
    if not completed_steps:
        # 没有进度，无法continue
        return ConfirmDecision(
            decision='cancel',
            payload={'error': 'no_progress', 'message': '无执行进度，请先选择执行选项'}
        )
    
    # 根据完成步骤确定下一步
    if "featurecounts" in completed_steps:
        return ConfirmDecision(decision='continue_analysis')
    elif "star" in completed_steps:
        return ConfirmDecision(decision='continue_featurecounts')  
    elif "fastp" in completed_steps:
        return ConfirmDecision(decision='continue_star')
    else:
        # 默认继续执行
        return ConfirmDecision(decision='continue')


def _handle_re_opt_command(
    completed_steps: List[str], 
    current_step: str, 
    raw_input: str
) -> ConfirmDecision:
    """处理重新优化命令"""
    
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
        decision=target_step,  # 直接返回目标步骤名
        execution_mode='optimized', 
        payload={'re_optimization': True, 'target_step': target_step, 'raw_input': raw_input}
    )


def get_available_commands(context: Dict[str, Any]) -> List[str]:
    """
    根据上下文获取当前可用的命令列表
    
    Args:
        context: 上下文信息
        
    Returns:
        可用命令列表
    """
    completed_steps = context.get('completed_steps', [])
    current_step = context.get('current_step', '')
    
    commands = []
    
    # 根据执行状态动态生成命令
    if completed_steps:
        commands.append('/continue')
        commands.append('/restart')
    else:
        commands.append('/execute_opt')
        commands.append('/yolo')  # 新增YOLO命令
    
    # 二次优化
    if current_step in {"fastp", "star", "featurecounts"}:
        commands.append('/re_opt')
    
    # 通用命令
    commands.extend(['/modify', '/cancel', '/quit'])
    
    return commands


def validate_command(command: str, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    验证命令是否在当前上下文中有效
    
    Args:
        command: 要验证的命令
        context: 上下文信息
        
    Returns:
        验证结果字典，包含 valid, reason 字段
    """
    available_commands = get_available_commands(context)
    
    # 检查命令是否在可用列表中
    normalized_cmd = command.lower().strip()
    
    # 检查直接匹配
    if normalized_cmd in [cmd.lower() for cmd in available_commands]:
        return {'valid': True, 'reason': 'direct_match'}
    
    # 检查别名匹配
    if normalized_cmd in COMMAND_ALIASES:
        return {'valid': True, 'reason': 'alias_match'}
    
    # 检查modify前缀
    for prefix in MODIFY_PREFIXES:
        if normalized_cmd.startswith(prefix):
            if '/modify' in [cmd.lower() for cmd in available_commands]:
                return {'valid': True, 'reason': 'modify_with_params'}
    
    return {
        'valid': False, 
        'reason': 'not_available',
        'available_commands': available_commands
    }