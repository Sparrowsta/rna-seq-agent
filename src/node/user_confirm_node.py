"""用户确认节点 - 纯数字选择模式

展示配置并等待用户纯数字选择决策。
完全弃用斜杠命令，采用数字索引选择模式。
"""

import os
from typing import Dict, Any
from ..state import AgentState
"""交互输出统一使用标准 print，避免额外封装"""
from .confirm import (
    build_confirm_view, render_confirm, 
    parse_numeric_selection, get_execution_mode_selection, parse_execution_mode_selection
)


async def user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """用户确认节点 - 展示配置并等待用户纯数字选择"""
    
    try:
        return await _numeric_user_confirm_node(state)
    except Exception as e:
        print(f"❌ 确认节点出错: {e}")
        # 返回安全的取消状态
        return {
            "response": f"确认节点执行错误: {e}",
            "status": "error",
            "user_decision": "cancel"
        }


async def _numeric_user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """纯数字选择用户确认节点实现"""
    
    # 1. 构建视图模型
    view = build_confirm_view(state)
    
    # 2. 渲染输出
    rendered_lines = render_confirm(view)
    for line in rendered_lines:
        print(line)
    
    # 3. 循环获取用户输入直到有效
    while True:
        try:
            user_choice = input("请输入数字选择 (0取消): ").strip()
            
            # 构建命令解析上下文
            context = {
                'completed_steps': getattr(state, 'completed_steps', []),
                'current_step': getattr(state, 'current_step', ''),
                'batch_optimizations': getattr(state, 'batch_optimizations', {}),
                'batch_optimization_complete': getattr(state, 'batch_optimization_complete', False)
            }
            
            # 使用新的数字选择解析
            decision = parse_numeric_selection(user_choice, view.commands, context)
            
            # 4. 检查是否有错误
            if decision.payload.get('error'):
                error_message = decision.payload.get('message', '输入错误')
                print(f"❌ {error_message}")
                continue  # 重新提示用户输入
            
            # 5. 处理需要模式选择的execute决策
            if decision.decision == 'execute' and decision.payload.get('needs_mode_selection'):
                execution_mode = await _handle_execution_mode_selection()
                if execution_mode == 'cancel':
                    continue  # 用户取消了，回到主菜单
                decision.execution_mode = execution_mode
            
            # 6. 处理需要二次输入的modify决策
            if decision.needs_modify_content:
                modify_content = await _handle_modify_content_input()
                if modify_content is None:
                    continue  # 用户取消了修改，回到主菜单
                decision.modify_content = modify_content
            
            # 7. 生成决策消息并返回结果
            decision_message = _generate_decision_message(decision)
            print(f"\n✅ {decision_message}")
            
            return _build_node_result(state, decision)
            
        except KeyboardInterrupt:
            print("\n\n❌ 用户中断，返回普通模式")
            return {
                "response": "用户中断操作",
                "status": "cancel",
                "user_decision": "cancel",
                "routing_decision": "normal"
            }
        except Exception as e:
            print(f"❌ 输入处理错误: {e}")
            continue


async def _handle_execution_mode_selection() -> str:
    """处理执行模式选择"""
    print("\n" + "\n".join(get_execution_mode_selection()))
    
    while True:
        try:
            mode_input = input("请选择执行模式: ").strip()
            execution_mode = parse_execution_mode_selection(mode_input)
            
            if execution_mode is None:
                print("❌ 无效选择，请输入 1-3 选择模式，或 0 返回")
                continue
            
            if execution_mode == 'cancel':
                return 'cancel'  # 用户选择返回
            
            # 显示选择的模式
            mode_descriptions = {
                'single': '单次执行 - 直接执行不优化',
                'optimized': '优化模式 - 每步优化后确认',
                'batch_optimize': '批次优化 - 收集所有优化建议后统一处理'
            }
            
            description = mode_descriptions.get(execution_mode, execution_mode)
            print(f"✅ 已选择: {description}")
            return execution_mode
            
        except KeyboardInterrupt:
            return 'cancel'
        except Exception as e:
            print(f"❌ 模式选择错误: {e}")
            continue


async def _handle_modify_content_input() -> str:
    """处理修改内容二次输入"""
    try:
        modify_input = input("请输入修改内容 (回车取消): ").strip()
        
        if not modify_input:
            print("❌ 修改取消")
            return None  # 空输入表示取消
        
        print(f"✅ 修改内容: {modify_input}")
        return modify_input
        
    except KeyboardInterrupt:
        print("\n❌ 修改取消")
        return None


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
        mode = decision.execution_mode or 'unknown'
        return f"⚡ 执行RNA-seq流水线 ({mode}模式)"
    elif decision.is_continue:
        return "➡️ 继续到下一步"
    elif decision.is_workflow_step:
        return f"🎯 执行 {decision.decision.upper()} 步骤"
    elif decision.decision == 'modify':
        content_info = f" - {decision.modify_content}" if decision.modify_content else ""
        return f"🔧 修改配置{content_info}"
    elif decision.decision == 'cancel':
        return "❌ 取消分析"
    elif decision.decision == 'quit':
        return "🚪 退出程序"
    elif decision.decision == 'restart':
        return "🔄 重新开始"
    else:
        return f"🎯 {decision.decision}"


def _build_node_result(state: AgentState, decision) -> Dict[str, Any]:
    """构建节点返回结果"""
    
    # 基础返回字段
    result = {
        "response": _generate_decision_message(decision),
        "user_decision": decision.decision,
        "status": "confirm_complete"
    }
    
    # 添加执行模式
    if decision.execution_mode:
        result["execution_mode"] = decision.execution_mode
    
    # 添加修改内容
    if decision.modify_content:
        result["modify_requirements"] = {
            "raw_input": decision.modify_content,
            "source": "numeric_selection"
        }
    
    # 添加payload中的特殊字段
    if decision.payload.get('restart'):
        result["restart_requested"] = True
    
    if decision.payload.get('re_optimization'):
        result["re_optimization_target"] = decision.payload.get('target_step')
    
    # 处理特殊路由
    if decision.decision == 'cancel':
        result["routing_decision"] = "normal"
    elif decision.decision == 'quit':
        result["routing_decision"] = "end"
    
    return result
