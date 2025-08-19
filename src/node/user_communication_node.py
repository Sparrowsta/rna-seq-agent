from typing import Dict, Any
from ..state import UserCommunicationNodeState

async def user_communication_node(state: UserCommunicationNodeState) -> Dict[str, Any]:
    """User Communication节点 - 用户交互入口"""
    print(f"🔬 RNA-seq智能分析助手")
    print(f"   输入 /plan 开始分析 | /help 获取帮助 | /exit 退出")
    
    # 检查并显示来自normal节点的结果
    if hasattr(state, 'query_response') and state.query_response:
        print()
        print(f"🎯 {state.query_response}")
        print()
    
    # 获取用户输入
    try:
        user_input = input("请输入: ").strip()
        print(f"📝 收到输入: {user_input}")
        
        # 基本路由判断
        if user_input.lower() in ['/exit', 'exit', '退出']:
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "end",
                "response": "再见！",
                "status": "exiting"
            }
        elif user_input.lower() in ['/plan', '开始分析']:
            return {
                "messages": [{"role": "user", "content": user_input}], 
                "routing_decision": "plan",
                "response": "进入分析计划流程...",
                "status": "routing_to_plan"
            }
        else:
            # 其他输入交给normal节点处理
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "normal", 
                "response": "正在分析您的需求...",
                "status": "routing_to_normal"
            }
            
    except KeyboardInterrupt:
        return {
            "messages": [{"role": "user", "content": "KeyboardInterrupt"}],
            "routing_decision": "end",
            "response": "用户中断退出",
            "status": "interrupted"
        }