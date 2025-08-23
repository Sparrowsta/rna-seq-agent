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
        elif user_choice_lower in ['/replan', '/重新规划', '/修改'] or user_choice_lower.startswith('/replan '):
            user_decision = "replan"
            decision_msg = "🔄 重新规划配置"
            
            # 处理/replan命令中的新需求
            if user_choice_lower.startswith('/replan'):
                # 智能提取/replan后面的内容，处理有无空格的情况
                replan_content = user_choice_lower.replace('/replan', '', 1).strip()
                if replan_content:
                    print(f"📝 检测到新配置需求: {replan_content}")
                    # 将新需求直接保存为user_requirements，让Plan节点的LLM来解析
                    new_user_requirements = {"raw_input": replan_content}
                else:
                    # 纯/replan命令，清空旧需求
                    new_user_requirements = {}
                    print(f"📝 清空旧配置需求，重新规划")
        elif user_choice_lower in ['/cancel', '/取消']:
            user_decision = "cancel"
            decision_msg = "❌ 取消分析"
        else:
            # 无效输入，提示用户重新选择
            print(f"❌ 无效输入: {user_choice}")
            print(f"请选择有效的命令: /execute, /replan, /cancel, /quit")
            # 递归调用自己，重新获取用户输入
            return await user_confirm_node(state)
        
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
        
        # 重新规划时设置replan需求，保持初始user_requirements不变
        "user_requirements": getattr(state, 'user_requirements', {}),  # 保持初始需求
        "replan_requirements": new_user_requirements if 'new_user_requirements' in locals() else {},  # replan需求
        
        # 保存用户选择用于后续处理
        "messages": [{"role": "user", "content": user_choice}]
    }