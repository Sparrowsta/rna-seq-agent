#!/usr/bin/env python3
"""
RNA-seq智能分析助手 - 主程序入口
基于LangGraph Plan-and-Execute架构的AI Agent系统
"""

import os
import sys
import asyncio
from typing import Dict, Any
from pathlib import Path
# 确保项目根目录在Python路径中
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from src.state import NormalNodeState
from src.graph import create_agent
from src.core import test_llm_connection

# 导入必要的组件
from dotenv import load_dotenv

def load_environment():
    """加载环境变量配置"""
    env_path = project_root / "config" / ".env"
    
    if env_path.exists():
        load_dotenv(env_path)
        print(f"✅ 已加载环境配置: {env_path}")
    else:
        print("⚠️  未找到环境配置文件: config/.env")
    
    # 验证必要的环境变量
    if not os.environ.get("DEEPSEEK_API_KEY"):
        print("❌ 错误: 未找到DEEPSEEK_API_KEY环境变量")
        print("请在config/.env文件中设置: DEEPSEEK_API_KEY=your-api-key")
        sys.exit(1)


def create_deepseek_llm():
    """创建并测试DeepSeek LLM实例"""
    success, message = test_llm_connection()
    
    if success:
        print(f"✅ DeepSeek LLM连接成功: {message}")
    else:
        print(f"❌ DeepSeek LLM连接失败: {message}")
        sys.exit(1)
        
    return success

async def run_interactive_session(agent):
    """运行交互式会话"""
    print("\n💬 RNA-seq智能分析助手启动")
    print("🔹 系统将直接进入用户通信模式")
    print("🔹 Agent将处理所有用户交互\n")
    
    # 创建空的初始状态，让user_communication节点来处理输入
    initial_state = {
        "response": "",
        "status": "starting"
    }
    
    try:
        # 调用Agent - 从user_communication节点开始
        result = await agent.ainvoke(initial_state)
        print("🤖 会话结束")
        
    except KeyboardInterrupt:
        print("\n👋 收到中断信号，退出程序")
    except Exception as e:
        print(f"❌ 处理错误: {e}")


async def main():
    """主函数"""
    # 加载环境配置
    load_environment()
    
    # 测试LLM连接
    create_deepseek_llm()
    
    # 创建Agent
    agent = create_agent()
    
    # 运行交互式会话
    await run_interactive_session(agent)


if __name__ == "__main__":
    asyncio.run(main())