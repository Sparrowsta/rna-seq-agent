"""
RNA-seq智能分析助手 - 核心包
基于LangGraph的AI Agent系统
"""

__version__ = "1.0.0"
__author__ = "RNA-seq Team"
__description__ = "基于AI的生物信息学RNA-seq分析助手"

# 导出核心组件
from .state import AgentState
from .graph import create_agent
from .core import get_shared_llm

__all__ = [
    "AgentState",
    "create_agent", 
    "get_shared_llm"
]