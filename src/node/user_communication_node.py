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
        
        # 基本路由判断
        if user_input.lower() in ['/exit', 'exit', '退出']:
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "end",
                "response": "再见！",
                "status": "normal"
            }
        elif user_input.lower().startswith('/plan'):
            # 智能提取/plan后面的内容，处理有无空格的情况
            plan_content = user_input.replace('/plan', '', 1).strip()
            
            if plan_content:
                print(f"📝 检测到分析需求: {plan_content}")
                # 将新需求直接保存为user_requirements，让Plan节点的LLM来解析
                plan_user_requirements = {"raw_input": plan_content}
                response_msg = f"进入分析计划流程...\n📝 分析需求: {plan_content}"
            else:
                # 纯/plan命令，无额外需求
                plan_user_requirements = {}
                response_msg = "进入分析计划流程..."
                print("📝 纯/plan命令，无额外需求")
            
            return {
                "messages": [{"role": "user", "content": user_input}], 
                "routing_decision": "plan",
                "response": response_msg,
                "user_requirements": plan_user_requirements,  # 传递给plan节点
                "status": "plan"
            }
        elif user_input.lower() in ['开始分析']:
            return {
                "messages": [{"role": "user", "content": user_input}], 
                "routing_decision": "plan",
                "response": "进入分析计划流程...",
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