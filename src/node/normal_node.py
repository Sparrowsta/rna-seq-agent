import os
from typing import Dict, Any
from ..state import NormalNodeState

from langchain_deepseek import ChatDeepSeek
from langchain_core.messages import HumanMessage
    
def create_normal_llm():
    """创建Normal节点专用的结构化LLM"""
    llm = ChatDeepSeek(
        model="deepseek-chat",
        api_key=os.environ["DEEPSEEK_API_KEY"],
        temperature=0.1
    )
    return llm.with_structured_output(NormalNodeState, method="json_mode")

async def normal_node(state: NormalNodeState) -> Dict[str, Any]:
    """Normal节点 - 用户交互和信息查询"""
    print(f"💬 处理用户请求...")
    print(f"   用户输入: {state.get('input', '')}")
    
    user_input = state.get("input", "").lower()
    
    # 构建LLM提示
    prompt = f"""
你是RNA-seq智能分析助手的Normal模式处理器。请分析用户输入并返回JSON格式响应。

用户输入: "{user_input}"

分析任务:
1. 识别查询类型 (info/analysis/help/plan)
2. 决定路由方向 (normal/plan) 
3. 解析用户意图
4. 生成查询响应和建议操作
5. 提供用户友好的回复

查询类型说明:
- info: 信息查询 (如查看文件、基因组列表、帮助信息)
- analysis: 分析意图识别 (用户想进行RNA-seq分析)
- help: 帮助请求 (如何使用、功能介绍)  
- plan: 明确进入计划模式请求 (如"/plan", "/开始分析")

路由决策逻辑:
- normal: 继续在Normal模式处理 (信息查询、帮助)
- plan: 进入Plan模式 (用户表达分析意图或明确要求)

请返回包含所有必需字段的JSON格式响应。
"""
    
    try:
        # 调用结构化LLM
        structured_llm = create_normal_llm()
        response = await structured_llm.ainvoke([HumanMessage(content=prompt)])
        
        print(f"🤖 LLM分析结果: {response.query_type}, 路由: {response.routing_decision}")
        
        # 将Pydantic模型转换为字典返回
        return response.dict()
        
    except Exception as e:
        print(f"❌ Normal节点LLM调用失败: {e}")
        return {
            "query_type": "help",
            "routing_decision": "normal",
            "query_response": f"处理请求时出现错误: {str(e)}",
            "user_intent": "系统错误",
            "suggested_actions": ["请重试或联系管理员"],
            "response": "抱歉，处理您的请求时遇到问题，请重试",
            "status": "error"
        }