"""
LLM管理模块 - 重构版
集中管理DeepSeek LLM实例，支持配置系统
"""

from typing import Optional
from langchain_deepseek import ChatDeepSeek
from langchain_core.messages import HumanMessage
from .config.settings import Settings
from .logging_bootstrap import get_logger, log_llm_preview

logger = get_logger("rna.core")

class LLMManager:
    """LLM管理器 - 负责创建和管理LLM实例"""
    
    def __init__(self, settings: Settings):
        self.settings = settings
        self._llm_instance: Optional[ChatDeepSeek] = None
    
    def get_llm(self) -> ChatDeepSeek:
        """获取LLM实例"""
        if self._llm_instance is None:
            self._llm_instance = ChatDeepSeek(
                model="deepseek-chat",
                api_key=self.settings.deepseek_api_key,
                temperature=self.settings.llm_temperature
            )
        return self._llm_instance
    
    def test_connection(self) -> tuple[bool, str]:
        """测试LLM连接"""
        try:
            llm = self.get_llm()
            test_response = llm.invoke([HumanMessage(content="测试连接，请回复'连接成功'")])
            # 仅记录返回内容的预览（不记录提示/指令）
            log_llm_preview(logger, "core.test", getattr(test_response, "content", test_response))
            # 处理不同的响应类型
            content = test_response.content
            if isinstance(content, str):
                return True, content
            elif isinstance(content, list):
                return True, str(content)
            else:
                return True, "连接成功"
        except Exception as e:
            return False, str(e)

# ==================== 向后兼容性支持 ====================
# 保留原有的全局函数接口，以免破坏现有代码

_global_llm_manager: Optional[LLMManager] = None

def _get_global_manager() -> LLMManager:
    """获取全局LLM管理器实例"""
    global _global_llm_manager
    if _global_llm_manager is None:
        # 创建默认配置
        settings = Settings()
        _global_llm_manager = LLMManager(settings)
    return _global_llm_manager

def get_shared_llm() -> ChatDeepSeek:
    """获取共享的DeepSeek LLM实例 - 向后兼容接口"""
    manager = _get_global_manager()
    return manager.get_llm()

def test_llm_connection() -> tuple[bool, str]:
    """测试LLM连接 - 向后兼容接口"""
    manager = _get_global_manager()
    return manager.test_connection()