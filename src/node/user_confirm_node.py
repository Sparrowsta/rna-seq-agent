"""用户确认节点 - 重构版本

展示配置并等待用户决策，采用模块化架构降低复杂度。
使用 confirm 子模块提供结构化的视图、渲染和命令解析。
"""

import os
from typing import Dict, Any
from ..state import AgentState
from .confirm import (
    build_confirm_view, render_confirm, parse_confirm_command,
    get_available_commands, validate_command
)


async def user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """用户确认节点 - 展示配置并等待用户决策"""
    
    # 检查是否启用新版本（支持回退）
    use_new_version = os.getenv('CONFIRM_V2', 'true').lower() == 'true'
    
    if not use_new_version:
        # 回退到原始实现
        return await _legacy_user_confirm_node(state)
    
    try:
        return await _new_user_confirm_node(state)
    except Exception as e:
        print(f"❌ 新版确认节点出错，回退到原版: {e}")
        return await _legacy_user_confirm_node(state)


def _generate_decision_message(decision) -> str:
    """生成决策消息"""
    # 检查是否有预设的模式描述
    if decision.payload.get('mode_description'):
        return decision.payload['mode_description']
    
    # 检查是否为重新优化命令
    if decision.payload.get('re_optimization'):
        target_step = decision.decision.upper()
        return f"♻️ 重新优化当前步骤: {target_step}"
    
    # 根据决策类型生成消息
    if decision.is_execute:
        return "⚡ 执行RNA-seq流水线"
    elif decision.is_continue:
        return "➡️ 继续到下一步"
    elif decision.is_workflow_step:
        return f"🎯 执行 {decision.decision.upper()} 步骤"
    elif decision.decision == 'modify':
        return "🔧 修改配置"
    elif decision.decision == 'cancel':
        return "❌ 取消分析"
    elif decision.decision == 'quit':
        return "🚪 退出程序"
    elif decision.decision == 'restart':
        return "🔄 重新开始"
    else:
        return f"🎯 {decision.decision}"


async def _new_user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """新版用户确认节点实现"""
    
    # 1. 构建视图模型
    view = build_confirm_view(state)
    
    # 2. 渲染输出
    rendered_lines = render_confirm(view)
    for line in rendered_lines:
        print(line)
    
    # 3. 获取用户输入并解析
    user_choice = ""
    try:
        user_choice = input("请输入命令: ").strip()
        
        # 构建命令解析上下文
        context = {
            'completed_steps': getattr(state, 'completed_steps', []),
            'current_step': getattr(state, 'current_step', ''),
            'batch_optimizations': getattr(state, 'batch_optimizations', {}),
            'batch_optimization_complete': getattr(state, 'batch_optimization_complete', False)
        }
        
        # 解析命令
        decision = parse_confirm_command(user_choice, context)
        
        # 4. 处理特殊决策（如execute需要模式选择）
        if decision.decision == 'execute' and decision.payload.get('needs_mode_selection'):
            decision = await _handle_execution_mode_selection(decision, user_choice)
        
        # 5. 处理无效命令
        if decision.decision == 'cancel' and decision.payload.get('error') == 'invalid_command':
            return await _handle_invalid_command(state, user_choice, context)
        
        # 6. 构建返回结果
        return _build_node_result(state, decision, view)
        
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断，取消分析")
        decision_msg = "❌ 用户中断取消"
        return _build_cancel_result(state, decision_msg, "/cancel")
    except Exception as e:
        print(f"❌ 输入处理错误: {e}")
        decision_msg = "❌ 输入错误取消"
        return _build_cancel_result(state, decision_msg, "/cancel")


async def _handle_execution_mode_selection(decision, original_input: str):
    """处理执行模式选择"""
    print(f"\n🔄 **请选择执行模式:**")
    print(f"   1. 单次执行 - 运行完整流水线，无优化处理")
    print(f"   2. 精细优化 - 运行流水线，每步立即应用优化参数") 
    print(f"   3. 批次优化 - 收集所有工具优化建议后统一处理")
    print(f"   0. 返回上级菜单")
    
    try:
        mode_choice = input("请选择执行模式 (1-3, 0返回): ").strip()
        
        if mode_choice == "1":
            decision.execution_mode = 'single'
            decision.payload['mode_description'] = "✅ 单次执行完整RNA-seq流水线"
        elif mode_choice == "2":
            decision.execution_mode = 'optimized'
            decision.payload['mode_description'] = "⚡ 精细优化执行完整RNA-seq流水线"
        elif mode_choice == "3":
            decision.execution_mode = 'batch_optimize'
            decision.payload['mode_description'] = "📦 批次优化执行完整RNA-seq流水线"
        elif mode_choice == "0":
            decision.decision = 'cancel'
            decision.payload['mode_description'] = "返回上级菜单"
        else:
            decision.decision = 'cancel'
            decision.payload['error'] = 'invalid_mode_selection'
            decision.payload['mode_description'] = f"❌ 无效选择: {mode_choice}"
            
    except KeyboardInterrupt:
        decision.decision = 'cancel'
        decision.payload['mode_description'] = "返回上级菜单"
    
    return decision


async def _handle_invalid_command(state: AgentState, user_choice: str, context: Dict[str, Any]):
    """处理无效命令"""
    print(f"❌ 无效输入: {user_choice}")
    
    # 获取并显示可用命令
    available_commands = get_available_commands(context)
    print(f"请选择有效的命令: {', '.join(available_commands)}")
    
    # 递归调用自己，重新获取用户输入
    return await user_confirm_node(state)


def _build_node_result(state: AgentState, decision, view) -> Dict[str, Any]:
    """构建节点返回结果"""
    
    # 获取基础配置
    nextflow_config = state.nextflow_config or {}
    resource_config = state.resource_config or {}
    config_reasoning = getattr(state, 'config_reasoning', '') or ''
    
    # 生成决策消息
    decision_msg = _generate_decision_message(decision)
    
    print(f"🎯 {decision_msg}")
    
    # 构建确认消息
    reasoning_line = f"💭 决策理由: {config_reasoning}\\n" if config_reasoning else ""
    confirmation_message = f"""🎯 分析配置已确认！

📋 配置项: {len(nextflow_config)} 个参数已设置
{reasoning_line}🎯 用户选择: {decision_msg}

准备进入下一阶段..."""

    # 构建返回字典
    result = {
        # 从 prepare_node 继承并传递
        "nextflow_config": nextflow_config,
        "resource_config": resource_config,
        "config_reasoning": config_reasoning,
        
        # 当前节点输出
        "confirmation_message": confirmation_message,
        "user_decision": decision.decision,
        "response": decision_msg,
        "status": decision.decision,
        "execution_mode": decision.execution_mode or getattr(state, 'execution_mode', 'single'),
        
        # 进度信息
        "completed_steps": getattr(state, 'completed_steps', []),
        "current_step": getattr(state, 'current_step', ''),
        
        # 批次优化相关状态
        "batch_optimizations": getattr(state, 'batch_optimizations', {}),
        "batch_optimization_complete": getattr(state, 'batch_optimization_complete', False),
        "batch_optimization_mode": (decision.execution_mode == 'batch_optimize'),
        "batch_optimization_round": getattr(state, 'batch_optimization_round', 1) + (1 if decision.decision == "execute" and decision.execution_mode == 'batch_optimize' else 0),
        
        # modify 需求处理
        "user_requirements": getattr(state, 'user_requirements', {}),
        "modify_requirements": {"raw_input": decision.modify_content or ""} if decision.modify_content else {},
        
        # 消息和Base快照持久化
        "messages": [{"role": "user", "content": decision.payload.get('raw_input', '')}],
        "prepare_defaults_nextflow_config": getattr(state, 'prepare_defaults_nextflow_config', {}),
        "prepare_defaults_resource_config": getattr(state, 'prepare_defaults_resource_config', {}),
        "prepare_defaults_fastp_params": getattr(state, 'prepare_defaults_fastp_params', {}),
        "prepare_defaults_star_params": getattr(state, 'prepare_defaults_star_params', {}),
        "prepare_defaults_featurecounts_params": getattr(state, 'prepare_defaults_featurecounts_params', {})
    }
    
    return result


def _build_cancel_result(state: AgentState, decision_msg: str, user_choice: str) -> Dict[str, Any]:
    """构建取消/错误结果"""
    nextflow_config = state.nextflow_config or {}
    config_reasoning = getattr(state, 'config_reasoning', '') or ''
    
    reasoning_line = f"💭 决策理由: {config_reasoning}\\n" if config_reasoning else ""
    confirmation_message = f"""🎯 分析配置已确认！

📋 配置项: {len(nextflow_config)} 个参数已设置
{reasoning_line}🎯 用户选择: {decision_msg}

准备进入下一阶段..."""
    
    return {
        "nextflow_config": nextflow_config,
        "resource_config": state.resource_config or {},
        "config_reasoning": config_reasoning,
        "confirmation_message": confirmation_message,
        "user_decision": "cancel",
        "response": decision_msg,
        "status": "cancel",
        "execution_mode": getattr(state, 'execution_mode', 'single'),
        "completed_steps": getattr(state, 'completed_steps', []),
        "current_step": getattr(state, 'current_step', ''),
        "batch_optimizations": getattr(state, 'batch_optimizations', {}),
        "batch_optimization_complete": getattr(state, 'batch_optimization_complete', False),
        "batch_optimization_mode": False,
        "batch_optimization_round": getattr(state, 'batch_optimization_round', 1),
        "user_requirements": getattr(state, 'user_requirements', {}),
        "modify_requirements": {},
        "messages": [{"role": "user", "content": user_choice}],
        "prepare_defaults_nextflow_config": getattr(state, 'prepare_defaults_nextflow_config', {}),
        "prepare_defaults_resource_config": getattr(state, 'prepare_defaults_resource_config', {}),
        "prepare_defaults_fastp_params": getattr(state, 'prepare_defaults_fastp_params', {}),
        "prepare_defaults_star_params": getattr(state, 'prepare_defaults_star_params', {}),
        "prepare_defaults_featurecounts_params": getattr(state, 'prepare_defaults_featurecounts_params', {})
    }


async def _legacy_user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """原始实现 - 用于回退"""
    # 这里导入并调用原始实现
    # 为了简化，我们直接抛出异常提示用户
    raise NotImplementedError("Legacy version not implemented in refactor. Please set CONFIRM_V2=true to use new version.")