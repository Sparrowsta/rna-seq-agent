"""用户确认节点 - 纯数字选择模式

展示配置并等待用户纯数字选择决策。
完全弃用斜杠命令，采用数字索引选择模式。
"""

from typing import Dict, Any
from ..state import AgentState
from ..config.default_tool_params import DEFAULT_FASTP_PARAMS, DEFAULT_STAR_PARAMS, DEFAULT_FEATURECOUNTS_PARAMS, DEFAULT_HISAT2_PARAMS
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

    # 在显示界面前，根据返回原因决定是否清空状态
    if state.return_reason == "completed":
        print(f"🎉 {state.execution_mode}模式任务完成，重置状态准备新任务...")
        _reset_state_for_execution_mode(state, state.execution_mode, preserve_base_config=True)
        print("✅ 状态重置完成")
    elif state.return_reason == "batch_collect":
        print(f"📋 Batch优化模式：保留收集的优化建议")
    elif state.return_reason == "step_confirm":
        print(f"🔄 Optimized模式步骤确认：保留当前进度")
    elif state.return_reason == "failed":
        print(f"⚠️  执行失败：保留错误信息供分析")

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
                'current_step': getattr(state, 'current_step', '')
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

    # 如果是执行决策，可以添加额外的日志
    if decision.decision == 'execute' and decision.execution_mode:
        print(f"✅ 用户确认执行{decision.execution_mode}模式")

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


def _reset_state_for_execution_mode(state: AgentState, mode: str, preserve_base_config: bool = True) -> None:
    """
    根据执行模式重置状态字段

    Args:
        state: AgentState实例
        mode: 执行模式 (single/optimized/batch_optimize/yolo)
        preserve_base_config: 是否保留基础配置（prepare_node生成的配置）
    """
    mode = (mode or '').lower()

    if mode in ('single', 'optimized'):
        # Single和Optimized模式：完全重置，除了基础配置
        _reset_tool_params(state)
        _reset_optimization_fields(state)
        _reset_user_modifications(state)
        _reset_execution_results(state)

    elif mode == 'batch_optimize':
        # Batch_optimize模式：只清空执行结果，保留优化累积状态
        _reset_execution_results(state)

    elif mode == 'yolo':
        # Yolo模式：当前与single模式相同（可以后续调整）
        _reset_tool_params(state)
        _reset_optimization_fields(state)
        _reset_user_modifications(state)
        _reset_execution_results(state)

    # 根据参数决定是否保留基础配置
    if not preserve_base_config:
        state.nextflow_config = {}
        state.resource_config = {}
        state.config_reasoning = ""


def _reset_tool_params(state: AgentState) -> None:
    """重置工具参数到默认值"""
    state.fastp_params = DEFAULT_FASTP_PARAMS.copy()
    state.star_params = DEFAULT_STAR_PARAMS.copy()
    state.hisat2_params = DEFAULT_HISAT2_PARAMS.copy()
    state.featurecounts_params = DEFAULT_FEATURECOUNTS_PARAMS.copy()


def _reset_optimization_fields(state: AgentState) -> None:
    """重置优化相关字段"""
    # 清空优化建议
    state.fastp_optimization_suggestions = ""
    state.star_optimization_suggestions = ""
    state.hisat2_optimization_suggestions = ""
    state.featurecounts_optimization_suggestions = ""

    # 清空优化参数变更
    state.fastp_optimization_params = {}
    state.star_optimization_params = {}
    state.hisat2_optimization_params = {}
    state.featurecounts_optimization_params = {}

    # 清空优化历史（保留最近的状态管理）
    state.fastp_optimization_history = []
    state.star_optimization_history = []
    state.hisat2_optimization_history = []
    state.featurecounts_optimization_history = []


def _reset_user_modifications(state: AgentState) -> None:
    """重置用户修改相关字段"""
    state.modification_history = []
    state.modify_requirements = {}


def _reset_execution_results(state: AgentState) -> None:
    """重置执行结果字段"""
    state.fastp_results = {}
    state.star_results = {}
    state.hisat2_results = {}
    state.featurecounts_results = {}

    # 重置分析结果
    state.overall_summary = ""
    state.key_findings = []
    state.sample_health_assessment = ""
    state.quality_metrics_analysis = ""
    state.optimization_recommendations = []
    state.risk_warnings = []
    state.next_steps = []
