#!/usr/bin/env python3
"""
RNA-seq智能分析助手 - 主程序入口
基于LangGraph Plan-and-Execute架构的AI Agent系统

重构版本 - 清理了导入路径，使用新的配置管理系统
"""

import sys
import asyncio
from pathlib import Path
from typing import Dict, Any

# 添加项目根目录到Python路径 - 修复导入问题
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

# 现在使用正确的导入路径
from src.config.settings import Settings
from src.state import AgentState
from src.graph import create_agent
from src.core import test_llm_connection

def initialize_application() -> Settings:
    """初始化应用程序配置"""
    print("🚀 初始化RNA-seq智能分析助手...")
    
    # 创建配置实例
    settings = Settings()
    
    # 验证环境配置
    is_valid, errors = settings.validate_environment()
    
    if not is_valid:
        print("❌ 环境配置验证失败:")
        for error in errors:
            print(f"   - {error}")
        sys.exit(1)
    
    # 显示配置信息
    print(f"✅ 环境类型: {'容器环境' if settings.is_container_environment else '本地开发环境'}")
    print(f"✅ 配置目录: {settings.config_dir}")
    print(f"✅ 数据目录: {settings.data_dir}")
    
    # 验证关键文件
    if settings.genomes_config_path.exists():
        print(f"✅ 基因组配置文件存在: {settings.genomes_config_path}")
    else:
        print(f"⚠️ 基因组配置文件不存在: {settings.genomes_config_path}")
    
    return settings

def validate_llm_connection() -> bool:
    """验证LLM连接"""
    print("🔗 验证DeepSeek LLM连接...")
    
    success, message = test_llm_connection()
    
    if success:
        print(f"✅ DeepSeek LLM连接成功: {message}")
        return True
    else:
        print(f"❌ DeepSeek LLM连接失败: {message}")
        return False

async def run_interactive_session(agent, settings: Settings):
    """运行交互式会话"""
    print("\n💬 RNA-seq智能分析助手启动")
    print("🔹 系统将直接进入用户通信模式")
    print("🔹 Agent将处理所有用户交互")
    print(f"🔹 工作目录: {settings.project_root}")
    print()
    
    # 创建初始状态
    initial_state = AgentState(status="normal")
    
    try:
        # 使用流式调用处理用户交互
        async for chunk in agent.astream(initial_state):
            # 处理节点更新
            for node_name, node_update in chunk.items():
                if node_update and isinstance(node_update, dict):
                    # 显示响应信息
                    response = node_update.get("response", "")
                    if response:
                        print(f"🔄 [{node_name}]: {response}")
                    
                    # 显示状态更新
                    status = node_update.get("status", "")
                    if status and status != "normal":
                        print(f"📊 状态更新: {status}")
        
        print("🤖 会话正常结束")
        
    except KeyboardInterrupt:
        print("\n👋 收到中断信号，正在安全退出...")
    except Exception as e:
        print(f"❌ 运行时错误: {e}")
        if settings.debug_mode:
            import traceback
            traceback.print_exc()

async def main():
    """主函数 - 应用程序入口点"""
    try:
        # 1. 初始化配置
        settings = initialize_application()
        
        # 2. 验证LLM连接  
        if not validate_llm_connection():
            sys.exit(1)
        
        # 3. 创建Agent
        print("⚙️ 创建LangGraph Agent...")
        agent = create_agent()
        print("✅ Agent创建成功")
        
        # 4. 运行交互式会话
        await run_interactive_session(agent, settings)
        
    except Exception as e:
        print(f"❌ 启动失败: {e}")
        sys.exit(1)

if __name__ == "__main__":
    asyncio.run(main())