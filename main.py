#!/usr/bin/env python3
"""
RNA-seq分析Agent主程序
遵循单一职责原则：专门处理程序启动和命令行交互
"""

import os
import sys
import logging
from pathlib import Path
from dotenv import load_dotenv

# 加载环境变量
load_dotenv(dotenv_path='config/.env')

# 导入agent模块
from agent.graph import agent_executor, print_graph_info, validate_graph_structure
from agent.state import create_initial_state
from agent.ui_manager import get_ui_manager

# ============================================================================
# 日志配置 - 遵循配置分离原则
# ============================================================================

def setup_logging(log_level: str = "INFO", log_file: str = None):
    """
    设置日志配置
    
    应用KISS原则：简单的日志配置
    """
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # 设置日志级别
    level = getattr(logging, log_level.upper(), logging.INFO)
    
    # 配置日志处理器
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        # 确保日志目录存在
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    # 配置根日志器
    logging.basicConfig(
        level=level,
        format=log_format,
        handlers=handlers,
        force=True  # 覆盖现有配置
    )
    
    # 设置第三方库的日志级别
    logging.getLogger("httpx").setLevel(logging.WARNING)
    logging.getLogger("openai").setLevel(logging.WARNING)
    logging.getLogger("langchain").setLevel(logging.WARNING)

def validate_system_requirements(silent=False):
    """
    验证系统要求
    
    应用KISS原则：简单的系统验证
    """
    if not silent:
        print("🔍 验证系统要求...")
    
    validation_results = []
    
    # 检查Python版本
    python_version = sys.version_info
    if python_version >= (3, 8):
        validation_results.append(("✅", f"Python版本: {python_version.major}.{python_version.minor}.{python_version.micro}"))
    else:
        validation_results.append(("❌", f"Python版本过低: {python_version.major}.{python_version.minor}.{python_version.micro} (需要 >= 3.8)"))
    
    # 检查必需的环境变量
    required_env_vars = ["DEEPSEEK_API_KEY"]
    for var in required_env_vars:
        if os.environ.get(var):
            validation_results.append(("✅", f"环境变量 {var}: 已设置"))
        else:
            validation_results.append(("❌", f"环境变量 {var}: 未设置"))
    
    # 检查配置文件
    config_files = [
        "config/genomes.json",
        "config/nextflow.config",
        "main.nf"
    ]
    
    for config_file in config_files:
        if os.path.exists(config_file):
            validation_results.append(("✅", f"配置文件 {config_file}: 存在"))
        else:
            validation_results.append(("❌", f"配置文件 {config_file}: 不存在"))
    
    # 检查图结构
    try:
        if validate_graph_structure():
            validation_results.append(("✅", "Agent图结构: 验证通过"))
        else:
            validation_results.append(("❌", "Agent图结构: 验证失败"))
    except Exception as e:
        validation_results.append(("❌", f"Agent图结构: 验证出错 - {str(e)}"))
    
    if not silent:
        # 显示验证结果
        print("\n📋 验证结果:")
        for status, message in validation_results:
            print(f"  {status} {message}")
        
        # 统计结果
        success_count = sum(1 for status, _ in validation_results if status == "✅")
        total_count = len(validation_results)
        
        print(f"\n📊 总结: {success_count}/{total_count} 项验证通过")
        
        if success_count == total_count:
            print("🎉 系统验证完全通过！")
            return True
        else:
            print("⚠️  系统验证存在问题，请检查上述失败项。")
            return False
    else:
        # 静默模式，只返回验证结果
        return validation_results

# ============================================================================
# 主程序入口 - 遵循命令模式
# ============================================================================

def main():
    """
    主程序入口
    
    应用模板方法模式：标准的程序启动流程
    """
    try:
        # 设置日志（默认WARNING级别）
        setup_logging("WARNING")
        
        logger = logging.getLogger(__name__)
        
        # 获取系统验证结果（静默模式）
        validation_results = validate_system_requirements(silent=True)
        
        # 获取UI管理器并显示欢迎信息
        ui_manager = get_ui_manager()
        ui_manager.show_welcome_banner(validation_results)
        
        # 创建初始状态
        initial_state = create_initial_state()
        
        # 直接运行agent
        try:
            final_state = agent_executor.invoke(initial_state, {"recursion_limit": 100})
        except KeyboardInterrupt:
            print("\n\n⚠️  用户中断程序")
            logger.info("用户中断程序")
        except Exception as e:
            print(f"\n❌ 程序执行出错: {str(e)}")
            logger.error(f"程序执行出错: {str(e)}")
    
    except KeyboardInterrupt:
        print("\n\n👋 程序被用户中断")
        sys.exit(0)
    
    except Exception as e:
        print(f"❌ 程序启动失败: {str(e)}")
        logging.error(f"程序启动失败: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()