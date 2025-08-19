from typing import Dict, Any
from ..state import UserConfirmState

async def user_confirm_node(state: UserConfirmState) -> Dict[str, Any]:
    """用户确认节点 - 展示配置并等待用户决策"""
    print(f"⏳ 等待用户确认...")
    
    # 展示当前配置摘要
    nextflow_config = state.get('nextflow_config', {})
    config_reasoning = state.get('config_reasoning', '')
    
    print(f"📋 配置摘要:")
    for key, value in nextflow_config.items():
        print(f"   {key}: {value}")
    
    print(f"\n💭 配置理由: {config_reasoning}")
    
    confirmation_message = f"""
🎯 分析配置已准备完成！
📋 配置项: {len(nextflow_config)} 个参数已设置
💭 决策理由: {config_reasoning}

请选择下一步操作:
- [E] Execute - 执行分析
- [M] Modify - 修改配置  
- [C] Cancel - 取消分析
"""
    
    # TODO: 实现真实的用户输入交互
    # 当前返回临时状态，等待后续集成交互逻辑
    
    return {
        # 从prepare_node继承并传递给execute_node
        "nextflow_config": nextflow_config,
        "config_reasoning": config_reasoning,
        
        # 当前节点输出
        "confirmation_message": confirmation_message,
        "user_decision": "",  # 等待用户输入
        "response": "配置已准备完成，等待用户确认",
        "status": "confirming"
    }