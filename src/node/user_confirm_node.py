from typing import Dict, Any
from ..state import AgentState

async def user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """用户确认节点 - 展示配置并等待用户决策"""
    print(f"\n{'='*60}")
    print(f"🎯 **分析配置确认**")
    print(f"{'='*60}")
    
    # 展示当前配置摘要
    nextflow_config = state.nextflow_config or {}
    config_reasoning = state.config_reasoning or "系统自动生成配置"
    
    print(f"\n📋 **配置摘要:**")
    if nextflow_config:
        for key, value in nextflow_config.items():
            # 格式化显示配置项
            if key == "genome_version":
                print(f"   🧬 基因组版本: {value}")
            elif key == "species":
                print(f"   🔬 物种: {value}")
            elif key == "qc_tool":
                print(f"   🧹 质控工具: {value}")
            elif key == "align_tool":
                print(f"   🎯 比对工具: {value}")
            elif key == "quant_tool":
                print(f"   📊 定量工具: {value}")
            else:
                print(f"   ⚙️ {key}: {value}")
    else:
        print(f"   ⚠️ 无配置信息")
    
    print(f"\n💭 **配置理由:**")
    print(f"   {config_reasoning}")
    
    print(f"\n🔄 **请选择下一步操作:**")
    print(f"   /execute  - 🚀 执行分析")
    print(f"   /replan   - 🔄 重新规划")  
    print(f"   /cancel   - ❌ 取消分析")
    print(f"   /quit     - 🚪 退出程序")
    print(f"{'='*60}")
    
    # 获取用户输入
    try:
        user_choice = input("请输入命令: ").strip()
        print(f"📝 用户输入: {user_choice}")
        
        # 处理用户输入 - 简化逻辑
        user_choice_lower = user_choice.lower()
        
        if user_choice_lower in ['/execute', '/执行']:
            user_decision = "execute"
            decision_msg = "✅ 确认执行分析"
        elif user_choice_lower in ['/quit', '/exit', 'quit', 'exit', '退出', 'bye']:
            user_decision = "quit"
            decision_msg = "🚪 退出程序"
        elif user_choice_lower in ['/replan', '/重新规划', '/修改']:
            user_decision = "replan"
            decision_msg = "🔄 重新规划配置"
        elif user_choice_lower in ['/cancel', '/取消']:
            user_decision = "cancel"
            decision_msg = "❌ 取消分析"
        else:
            # 其他所有输入都视为重新规划请求
            user_decision = "replan"
            decision_msg = "🔄 重新规划请求"
        
        print(f"🎯 {decision_msg}")
        
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断，取消分析")
        user_decision = "cancel"
        decision_msg = "❌ 用户中断取消"
    except Exception as e:
        print(f"❌ 输入处理错误: {e}")
        user_decision = "cancel"
        decision_msg = "❌ 输入错误取消"
    
    confirmation_message = f"""🎯 分析配置已确认！

📋 配置项: {len(nextflow_config)} 个参数已设置
💭 决策理由: {config_reasoning}
🎯 用户选择: {decision_msg}

准备进入下一阶段..."""
    
    return {
        # 从prepare_node继承并传递
        "nextflow_config": nextflow_config,
        "config_reasoning": config_reasoning,
        
        # 当前节点输出
        "confirmation_message": confirmation_message,
        "user_decision": user_decision,
        "response": decision_msg,
        "status": "confirm",
        
        # 将用户输入添加到messages中，供modify_node使用
        "messages": [{"role": "user", "content": user_choice}]
    }