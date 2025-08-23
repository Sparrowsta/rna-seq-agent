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

# 确保模块路径在Python路径中
sys.path.insert(0, "/src")

from src.state import AgentState
from src.graph import create_agent
from src.core import test_llm_connection

# 导入必要的组件

def load_environment():
    """加载环境变量配置并验证必要文件"""
    # 检查关键配置文件是否存在
    genomes_file = Path("/config/genomes.json")
    if genomes_file.exists():
        print(f"✅ 基因组配置文件存在: {genomes_file}")
    else:
        print(f"⚠️ 基因组配置文件不存在: {genomes_file}")
    
    # 检查Nextflow配置（如果有）
    nextflow_config = Path("/config/nextflow.config")
    if nextflow_config.exists():
        print(f"✅ Nextflow配置文件存在: {nextflow_config}")
    else:
        print(f"💡 Nextflow配置文件不存在: {nextflow_config} (可选)")
    
    # 验证环境变量（Docker --env-file 注入）
    if not os.environ.get("DEEPSEEK_API_KEY"):
        print("❌ 错误: 未找到DEEPSEEK_API_KEY环境变量")
        print("请确保config/.env文件存在且包含: DEEPSEEK_API_KEY=your-api-key")
        sys.exit(1)
    else:
        print("✅ 环境变量配置正确")


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
    
    # 创建完整的初始状态
    initial_state = AgentState(
        response="",
        status="normal"
    )
    
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