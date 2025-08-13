#!/usr/bin/env python3
"""
RNA-seq分析Agent主程序
遵循单一职责原则：专门处理程序启动和命令行交互
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from dotenv import load_dotenv

# 加载环境变量
load_dotenv(dotenv_path='config/.env')

# 导入agent模块
from agent.graph import agent_executor, print_graph_info, validate_graph_structure
from agent.state import create_initial_state
from agent.nodes.normal_mode_node import create_welcome_message

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

# ============================================================================
# 命令行参数解析 - 遵循命令模式
# ============================================================================

def create_argument_parser():
    """
    创建命令行参数解析器
    
    应用建造者模式：分步构建参数解析器
    """
    parser = argparse.ArgumentParser(
        description="RNA-seq分析Agent - 智能RNA测序数据分析助手",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python main.py                    # 启动交互模式
  python main.py --debug           # 启动调试模式
  python main.py --log-file logs/agent.log  # 指定日志文件
  python main.py --validate        # 验证系统配置
  python main.py --info            # 显示系统信息
        """
    )
    
    # 基本选项
    parser.add_argument(
        "--debug", 
        action="store_true",
        help="启用调试模式，显示详细日志"
    )
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="设置日志级别 (默认: INFO)"
    )
    
    parser.add_argument(
        "--log-file",
        type=str,
        help="指定日志文件路径"
    )
    
    # 系统选项
    parser.add_argument(
        "--validate",
        action="store_true",
        help="验证系统配置和依赖"
    )
    
    parser.add_argument(
        "--info",
        action="store_true",
        help="显示系统信息"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="RNA-seq Agent v1.0.0"
    )
    
    return parser

# ============================================================================
# 系统验证 - 遵循验证模式
# ============================================================================

def validate_system_requirements():
    """
    验证系统要求
    
    应用KISS原则：简单的系统验证
    """
    print("🔍 验证系统要求...")
    
    validation_results = []
    
    # 检查Python版本
    python_version = sys.version_info
    if python_version >= (3, 8):
        validation_results.append(("✅", f"Python版本: {python_version.major}.{python_version.minor}.{python_version.micro}"))
    else:
        validation_results.append(("❌", f"Python版本过低: {python_version.major}.{python_version.minor}.{python_version.micro} (需要 >= 3.8)"))
    
    # 检查必需的环境变量
    required_env_vars = ["OPENAI_API_KEY"]
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

def show_system_info():
    """
    显示系统信息
    
    应用信息专家模式：集中显示系统信息
    """
    print("=" * 60)
    print("🧬 RNA-seq分析Agent 系统信息")
    print("=" * 60)
    
    # 基本信息
    print(f"📍 工作目录: {os.getcwd()}")
    print(f"🐍 Python版本: {sys.version}")
    print(f"💻 操作系统: {os.name}")
    
    # 环境变量
    print("\n🔧 环境配置:")
    env_vars = ["OPENAI_API_KEY", "OPENAI_MODEL_NAME", "OPENAI_API_BASE"]
    for var in env_vars:
        value = os.environ.get(var, "未设置")
        if var == "OPENAI_API_KEY" and value != "未设置":
            value = f"{value[:8]}..." if len(value) > 8 else value
        print(f"  {var}: {value}")
    
    # 文件结构
    print("\n📁 项目结构:")
    important_paths = [
        "agent/",
        "config/",
        "main.nf",
        "main.py"
    ]
    
    for path in important_paths:
        if os.path.exists(path):
            if os.path.isdir(path):
                file_count = len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
                print(f"  📁 {path} ({file_count} 个文件)")
            else:
                size = os.path.getsize(path)
                print(f"  📄 {path} ({size} bytes)")
        else:
            print(f"  ❌ {path} (不存在)")
    
    print("=" * 60)

# ============================================================================
# 交互界面 - 遵循MVC模式
# ============================================================================

class InteractiveInterface:
    """
    交互界面控制器
    
    遵循单一职责原则：专门处理用户交互
    """
    
    def __init__(self, debug_mode: bool = False):
        self.debug_mode = debug_mode
        self.logger = logging.getLogger(__name__)
    
    def show_welcome(self):
        """显示欢迎信息"""
        print("\n" + "=" * 60)
        print("🧬 欢迎使用 RNA-seq分析Agent!")
        print("=" * 60)
        print("这是一个智能的RNA测序数据分析助手，可以帮助您：")
        print("• 📁 管理和查询FASTQ文件")
        print("• 🧬 配置基因组参考文件")
        print("• 📋 制定个性化分析计划")
        print("• 🚀 执行完整的RNA-seq分析流程")
        print("• 📊 生成分析结果报告")
        print("\n💡 提示：输入 'exit' 或 'quit' 退出程序")
        print("=" * 60)
    
    def show_goodbye(self):
        """显示告别信息"""
        print("\n" + "=" * 60)
        print("👋 感谢使用 RNA-seq分析Agent!")
        print("=" * 60)
        print("如果您有任何问题或建议，请联系技术支持。")
        print("祝您的研究工作顺利！🧬✨")
        print("=" * 60)
    
    def run_interactive_session(self):
        """
        运行交互会话
        
        应用状态机模式：管理会话状态
        """
        try:
            # 显示欢迎信息
            self.show_welcome()
            
            # 创建初始状态
            initial_state = create_initial_state()
            
            # 添加欢迎消息
            welcome_msg = create_welcome_message()
            initial_state["messages"] = [welcome_msg]
            
            if self.debug_mode:
                print_graph_info()
            
            # 运行agent
            self.logger.info("启动交互会话")
            
            try:
                # 增加递归限制配置以避免工具调用过多
                final_state = agent_executor.invoke(initial_state, {"recursion_limit": 100})
                
                if self.debug_mode:
                    print(f"\n[DEBUG] 最终状态: {final_state}")
                
            except KeyboardInterrupt:
                print("\n\n⚠️  用户中断程序")
                self.logger.info("用户中断程序")
            
            except Exception as e:
                print(f"\n❌ 程序执行出错: {str(e)}")
                self.logger.error(f"程序执行出错: {str(e)}")
                
                if self.debug_mode:
                    import traceback
                    traceback.print_exc()
            
            # 显示告别信息
            self.show_goodbye()
        
        except Exception as e:
            print(f"❌ 交互会话启动失败: {str(e)}")
            self.logger.error(f"交互会话启动失败: {str(e)}")
            sys.exit(1)

# ============================================================================
# 主程序入口 - 遵循命令模式
# ============================================================================

def main():
    """
    主程序入口
    
    应用模板方法模式：标准的程序启动流程
    """
    try:
        # 解析命令行参数
        parser = create_argument_parser()
        args = parser.parse_args()
        
        # 设置日志
        log_level = "DEBUG" if args.debug else args.log_level
        setup_logging(log_level, args.log_file)
        
        logger = logging.getLogger(__name__)
        logger.info("RNA-seq Agent 启动")
        
        # 处理特殊命令
        if args.info:
            show_system_info()
            return
        
        if args.validate:
            success = validate_system_requirements()
            sys.exit(0 if success else 1)
        
        # 验证系统要求（简化版）
        if not validate_system_requirements():
            print("\n⚠️  系统验证失败，但程序将继续运行。")
            print("某些功能可能无法正常工作。")
        
        # 启动交互界面
        interface = InteractiveInterface(debug_mode=args.debug)
        interface.run_interactive_session()
    
    except KeyboardInterrupt:
        print("\n\n👋 程序被用户中断")
        sys.exit(0)
    
    except Exception as e:
        print(f"❌ 程序启动失败: {str(e)}")
        logging.error(f"程序启动失败: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()