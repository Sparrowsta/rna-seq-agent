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

# 导入必要的组件
from dotenv import load_dotenv
from langchain_deepseek import ChatDeepSeek
from langchain_core.messages import HumanMessage

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
    """创建DeepSeek LLM实例"""

    llm = ChatDeepSeek(
            model="deepseek-chat",
            api_key=os.environ["DEEPSEEK_API_KEY"],
            temperature=0.1
        )
    
    try:
        # 测试连接
        test_response = llm.invoke([HumanMessage(content="测试连接，请回复'连接成功'")])
        print(f"✅ DeepSeek LLM连接成功: {test_response.content}")
        
        return llm
    
    except Exception as e:
        print(f"❌ DeepSeek LLM连接失败: {e}")
        sys.exit(1)

async def run_interactive_session(agent):
    """运行交互式会话"""
    print("\n💬 进入交互模式 (输入 'quit' 退出)")
    print("🔹 支持的命令:")
    print("   - 任何自然语言描述的RNA-seq分析需求")
    print("   - /plan - 制定分析计划")
    print("   - /execute - 执行分析")
    print("   - quit - 退出程序\n")
    
    state = NormalNodeState()
    
    while True:
        try:
            user_input = input("👤 您: ").strip()
            
            if user_input.lower() in ['quit', 'exit', '退出']:
                print("👋 再见！")
                break
            
            if not user_input:
                continue
            
            # 更新状态
            state["input"] = user_input
            state["messages"].append({"role": "user", "content": user_input})
            
            print("🤖 AI助手思考中...")
            
            # 调用Agent
            result = await agent.ainvoke(state)
            
            # 显示结果
            response = result.get("response", "处理完成")
            print(f"🤖 AI助手: {response}")
            
            # 更新状态
            state.update(result)
            
        except KeyboardInterrupt:
            print("\n👋 收到中断信号，退出程序")
            break
        except Exception as e:
            print(f"❌ 处理错误: {e}")


async def main():
    """主函数"""
    # 加载环境配置
    load_environment()
    
    # 创建LLM实例
    llm = create_deepseek_llm()
    
    # 创建Agent
    agent = create_agent()
    
    # 运行交互式会话
    await run_interactive_session(agent)


if __name__ == "__main__":
    asyncio.run(main())