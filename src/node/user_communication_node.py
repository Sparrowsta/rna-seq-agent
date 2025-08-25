from typing import Dict, Any
from ..state import AgentState

async def user_communication_node(state: AgentState) -> Dict[str, Any]:
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
        
        # 定义plan等价命令
        plan_prefixes = ['/plan', '/开始分析']
        user_input_lower = user_input.lower()
        is_plan_command = (user_input_lower in plan_prefixes or 
                          any(user_input_lower.startswith(f"{prefix} ") for prefix in plan_prefixes))
        
        # 基本路由判断
        if user_input_lower in ['/exit', '/退出']:
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "end",
                "response": "再见！",
                "status": "normal"
            }
        elif is_plan_command:
            # 优雅的参数提取 - 处理所有plan等价命令
            plan_content = ""
            for prefix in plan_prefixes:
                if user_input_lower.startswith(prefix.lower()):
                    plan_content = user_input[len(prefix):].strip()
                    break
            
            if plan_content:
                print(f"📝 检测到分析需求: {plan_content}")
                plan_user_requirements = {"raw_input": plan_content}
                response_msg = f"进入分析计划流程...\n📝 分析需求: {plan_content}"
            else:
                print("📝 纯plan命令，无额外需求")
                plan_user_requirements = {}
                response_msg = "进入分析计划流程..."
            
            return {
                "messages": [{"role": "user", "content": user_input}], 
                "routing_decision": "plan",
                "response": response_msg,
                "user_requirements": plan_user_requirements,
                "status": "plan"
            }
        else:
            # 其他输入交给normal节点处理
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "normal", 
                "response": "正在分析您的需求...",
                "status": "normal"
            }
            
    except KeyboardInterrupt:
        return {
            "messages": [{"role": "user", "content": "KeyboardInterrupt"}],
            "routing_decision": "end",
            "response": "用户中断退出",
            "status": "interrupted"
        }