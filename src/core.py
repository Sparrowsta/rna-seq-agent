"""
共享LLM实例模块
确保整个应用只使用一个DeepSeek LLM实例
"""
import os
from langchain_deepseek import ChatDeepSeek
from langchain_core.messages import HumanMessage

# 全局LLM实例
_llm_instance = None

def get_shared_llm():
    """获取共享的DeepSeek LLM实例"""
    global _llm_instance
    
    if _llm_instance is None:
        _llm_instance = ChatDeepSeek(
            model="deepseek-chat",
            api_key=os.environ["DEEPSEEK_API_KEY"],
            temperature=0.1
        )
        
    return _llm_instance

def test_llm_connection():
    """测试LLM连接"""
    try:
        llm = get_shared_llm()
        test_response = llm.invoke([HumanMessage(content="测试连接，请回复'连接成功'")])
        return True, test_response.content
    except Exception as e:
        return False, str(e)