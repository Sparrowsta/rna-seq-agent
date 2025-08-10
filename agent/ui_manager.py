"""
用户界面管理器
遵循单一职责原则：专门处理用户界面交互和UTF-8编码问题
"""

import sys
import os
import locale
import logging
from typing import Optional, Dict, Any, List
from datetime import datetime

# 配置日志
logger = logging.getLogger(__name__)

class EnhancedUIManager:
    """
    增强的用户界面管理器
    
    解决UTF-8编码问题并美化所有用户交互界面
    遵循单例模式：全局统一的UI管理
    """
    
    _instance = None
    _initialized = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        if not self._initialized:
            self._setup_encoding()
            self._setup_ui_components()
            self._initialized = True
    
    def _setup_encoding(self):
        """
        设置UTF-8编码环境
        
        应用KISS原则：简单有效的编码设置
        """
        try:
            # Windows系统设置
            if sys.platform.startswith('win'):
                os.system('chcp 65001 > nul 2>&1')
            
            # 设置locale
            try:
                locale.setlocale(locale.LC_ALL, 'zh_CN.UTF-8')
            except locale.Error:
                try:
                    locale.setlocale(locale.LC_ALL, 'C.UTF-8')
                except locale.Error:
                    logger.warning("无法设置UTF-8 locale")
            
            # 确保stdout/stderr使用UTF-8
            if hasattr(sys.stdout, 'reconfigure'):
                sys.stdout.reconfigure(encoding='utf-8', errors='ignore')
                sys.stderr.reconfigure(encoding='utf-8', errors='ignore')
            
            # 设置环境变量
            os.environ['PYTHONIOENCODING'] = 'utf-8'
            
        except Exception as e:
            logger.error(f"编码设置失败: {e}")
    
    def _setup_ui_components(self):
        """
        设置UI组件
        
        遵循策略模式：根据环境选择不同的UI策略
        """
        # 检查是否支持rich库
        self.use_rich = self._check_rich_support()
        
        if self.use_rich:
            self._setup_rich_ui()
        else:
            self._setup_basic_ui()
        
        # 设置颜色主题
        self.colors = {
            'primary': '\033[96m',      # 青色
            'secondary': '\033[92m',    # 绿色
            'warning': '\033[93m',      # 黄色
            'error': '\033[91m',        # 红色
            'info': '\033[94m',         # 蓝色
            'reset': '\033[0m',         # 重置
            'bold': '\033[1m',          # 粗体
            'dim': '\033[2m'            # 暗淡
        }
    
    def _check_rich_support(self) -> bool:
        """检查是否支持rich库"""
        try:
            import rich
            return True
        except ImportError:
            return False
    
    def _setup_rich_ui(self):
        """设置Rich UI组件"""
        try:
            from rich.console import Console
            from rich.prompt import Prompt
            from rich.panel import Panel
            from rich.text import Text
            from rich.table import Table
            from rich.progress import Progress
            
            self.console = Console(
                force_terminal=True,
                legacy_windows=False,
                width=100,
                color_system="auto"
            )
            
            self.rich_components = {
                'Prompt': Prompt,
                'Panel': Panel,
                'Text': Text,
                'Table': Table,
                'Progress': Progress
            }
            
        except ImportError as e:
            logger.warning(f"Rich库导入失败，使用基础UI: {e}")
            self.use_rich = False
            self._setup_basic_ui()
    
    def _setup_basic_ui(self):
        """设置基础UI组件"""
        self.console = None
        self.rich_components = {}
    
    def safe_print(self, text: str, style: str = 'reset', end: str = '\n'):
        """
        安全的打印函数，处理UTF-8编码问题
        
        应用装饰器模式：为打印功能添加编码处理
        """
        try:
            # 清理文本
            clean_text = self._clean_text(text)
            
            if self.use_rich and self.console:
                # 使用Rich打印
                if style in ['primary', 'secondary', 'warning', 'error', 'info']:
                    self.console.print(f"[{style}]{clean_text}[/]", end=end)
                else:
                    self.console.print(clean_text, end=end)
            else:
                # 使用基础打印
                color_code = self.colors.get(style, self.colors['reset'])
                reset_code = self.colors['reset']
                print(f"{color_code}{clean_text}{reset_code}", end=end, flush=True)
                
        except UnicodeEncodeError:
            # 编码失败时的回退方案
            fallback_text = text.encode('ascii', errors='ignore').decode('ascii')
            print(f"[编码问题] {fallback_text}", end=end, flush=True)
        except Exception as e:
            print(f"[打印错误] {str(e)}", end=end, flush=True)
    
    def _clean_text(self, text: str) -> str:
        """
        清理文本，处理编码问题
        
        应用KISS原则：简单有效的文本清理
        """
        try:
            if isinstance(text, bytes):
                text = text.decode('utf-8', errors='ignore')
            
            # 移除无效字符
            cleaned = text.encode('utf-8', errors='ignore').decode('utf-8')
            
            # 移除控制字符但保留换行符和制表符
            import re
            cleaned = re.sub(r'[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\x9F]', '', cleaned)
            
            return cleaned
        except Exception:
            return "文本包含无效字符，已清理"
    
    def show_welcome_banner(self):
        """
        显示欢迎横幅
        
        应用模板方法模式：标准的欢迎界面格式
        """
        banner = """
╔══════════════════════════════════════════════════════════════╗
║                    🧬 RNA-seq 分析助手                        ║
║                                                              ║
║              专业的RNA测序数据分析智能助手                      ║
║                                                              ║
║  功能特点：                                                   ║
║  • 智能对话式交互                                             ║
║  • 自动化流程管理                                             ║
║  • 多种基因组支持                                             ║
║  • Docker容器化部署                                           ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
        """
        
        self.safe_print(banner, 'primary')
        self.safe_print(f"启动时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 'dim')
        self.safe_print("=" * 60, 'dim')
        self.safe_print("")
    
    def get_user_input(self, prompt_text: str = "请输入您的需求", mode: str = "normal") -> str:
        """
        获取用户输入，处理UTF-8编码问题
        
        遵循策略模式：根据不同模式使用不同的提示样式
        """
        try:
            # 根据模式设置不同的提示符
            mode_icons = {
                'normal': '💬',
                'plan': '📋',
                'execute': '🚀'
            }
            
            icon = mode_icons.get(mode, '💬')
            
            if self.use_rich and self.console:
                # 使用Rich的输入
                from rich.prompt import Prompt
                user_input = Prompt.ask(
                    f"[cyan]{icon} {prompt_text}[/]",
                    console=self.console
                )
            else:
                # 使用基础输入
                self.safe_print(f"{icon} {prompt_text}: ", 'primary', end='')
                user_input = input()
            
            # 确保返回UTF-8字符串
            if isinstance(user_input, bytes):
                user_input = user_input.decode('utf-8', errors='ignore')
            
            return user_input.strip()
            
        except (UnicodeDecodeError, UnicodeEncodeError) as e:
            self.show_error(f"编码错误: {e}")
            return self.get_user_input(prompt_text, mode)
        
        except KeyboardInterrupt:
            self.show_info("用户取消操作")
            return "exit"
        
        except EOFError:
            self.show_info("输入结束")
            return "exit"
        
        except Exception as e:
            self.show_error(f"输入错误: {e}")
            return self.get_user_input(prompt_text, mode)
    
    def show_ai_response(self, response: str, mode: str = "normal"):
        """
        显示AI响应
        
        遵循装饰器模式：为AI响应添加美化装饰
        """
        # 清理响应文本
        clean_response = self._clean_text(response)
        
        if self.use_rich and self.console:
            # 使用Rich面板
            from rich.panel import Panel
            
            mode_titles = {
                'normal': '🤖 AI助手 - 信息收集模式',
                'plan': '📋 AI助手 - 计划制定模式', 
                'execute': '🚀 AI助手 - 执行模式'
            }
            
            title = mode_titles.get(mode, '🤖 AI助手')
            
            panel = Panel(
                clean_response,
                title=title,
                border_style="green",
                padding=(1, 2)
            )
            
            self.console.print(panel)
            self.console.print()
        else:
            # 使用基础格式
            self.safe_print("=" * 60, 'secondary')
            self.safe_print(f"🤖 AI助手 ({mode}模式)", 'secondary')
            self.safe_print("-" * 60, 'dim')
            self.safe_print(clean_response)
            self.safe_print("=" * 60, 'secondary')
            self.safe_print("")
    
    def show_error(self, message: str):
        """显示错误信息"""
        self.safe_print(f"❌ 错误: {message}", 'error')
    
    def show_warning(self, message: str):
        """显示警告信息"""
        self.safe_print(f"⚠️  警告: {message}", 'warning')
    
    def show_info(self, message: str):
        """显示信息"""
        self.safe_print(f"ℹ️  信息: {message}", 'info')
    
    def show_success(self, message: str):
        """显示成功信息"""
        self.safe_print(f"✅ 成功: {message}", 'secondary')
    
    def show_mode_switch(self, from_mode: str, to_mode: str, reason: str = ""):
        """
        显示模式切换信息
        
        应用观察者模式：通知用户模式变化
        """
        mode_names = {
            'normal': '信息收集模式',
            'plan': '计划制定模式',
            'execute': '执行模式'
        }
        
        from_name = mode_names.get(from_mode, from_mode)
        to_name = mode_names.get(to_mode, to_mode)
        
        switch_message = f"🔄 模式切换: {from_name} → {to_name}"
        if reason:
            switch_message += f"\n   原因: {reason}"
        
        if self.use_rich and self.console:
            from rich.panel import Panel
            panel = Panel(
                switch_message,
                title="模式切换",
                border_style="yellow",
                padding=(0, 2)
            )
            self.console.print(panel)
        else:
            self.safe_print("=" * 40, 'warning')
            self.safe_print(switch_message, 'warning')
            self.safe_print("=" * 40, 'warning')
        
        self.safe_print("")
    
    def show_progress_table(self, data: Dict[str, Any], title: str = "进度状态"):
        """
        显示进度表格
        
        应用表格模式：结构化显示进度信息
        """
        if self.use_rich and self.console:
            from rich.table import Table
            
            table = Table(title=title)
            table.add_column("步骤", style="cyan", width=20)
            table.add_column("状态", style="green", width=10)
            table.add_column("描述", style="dim", width=40)
            
            for step, info in data.items():
                if isinstance(info, dict):
                    status = "✅" if info.get("completed") else "⏳"
                    description = info.get("description", "")
                else:
                    status = "✅" if info else "⏳"
                    description = str(info) if info else ""
                
                table.add_row(step, status, description)
            
            self.console.print(table)
        else:
            # 基础表格格式
            self.safe_print(f"\n📊 {title}", 'info')
            self.safe_print("-" * 60, 'dim')
            
            for step, info in data.items():
                if isinstance(info, dict):
                    status = "✅" if info.get("completed") else "⏳"
                    description = info.get("description", "")
                else:
                    status = "✅" if info else "⏳"
                    description = str(info) if info else ""
                
                self.safe_print(f"{status} {step:<20} {description}", 'info')
            
            self.safe_print("-" * 60, 'dim')
        
        self.safe_print("")
    
    def confirm_action(self, message: str) -> bool:
        """
        确认操作
        
        应用KISS原则：简单的确认对话
        """
        try:
            if self.use_rich and self.console:
                from rich.prompt import Confirm
                return Confirm.ask(f"[yellow]{message}[/]", console=self.console)
            else:
                self.safe_print(f"❓ {message} (y/n): ", 'warning', end='')
                response = input().strip().lower()
                return response in ['y', 'yes', '是', '确认']
        
        except (KeyboardInterrupt, EOFError):
            return False
        except Exception as e:
            self.show_error(f"确认操作失败: {e}")
            return False
    
    def show_tool_call(self, tool_name: str, description: str = ""):
        """
        显示工具调用信息
        
        应用通知模式：通知用户工具调用状态
        """
        message = f"🔧 调用工具: {tool_name}"
        if description:
            message += f"\n   描述: {description}"
        
        self.safe_print(message, 'info')
    
    def show_separator(self, char: str = "=", length: int = 60, style: str = 'dim'):
        """显示分隔符"""
        self.safe_print(char * length, style)
    
    def clear_screen(self):
        """清屏"""
        try:
            os.system('cls' if os.name == 'nt' else 'clear')
        except Exception:
            # 如果清屏失败，打印空行
            self.safe_print("\n" * 50)

# ============================================================================
# 全局UI管理器实例 - 应用单例模式
# ============================================================================

_global_ui_manager: Optional[EnhancedUIManager] = None

def get_ui_manager() -> EnhancedUIManager:
    """
    获取全局UI管理器实例
    
    应用单例模式：确保全局只有一个UI管理器实例
    """
    global _global_ui_manager
    if _global_ui_manager is None:
        _global_ui_manager = EnhancedUIManager()
    return _global_ui_manager

# ============================================================================
# 便捷函数 - 遵循DRY原则
# ============================================================================

def safe_print(text: str, style: str = 'reset', end: str = '\n'):
    """全局安全打印函数"""
    get_ui_manager().safe_print(text, style, end)

def show_error(message: str):
    """全局错误显示函数"""
    get_ui_manager().show_error(message)

def show_success(message: str):
    """全局成功显示函数"""
    get_ui_manager().show_success(message)

def show_info(message: str):
    """全局信息显示函数"""
    get_ui_manager().show_info(message)

def show_warning(message: str):
    """全局警告显示函数"""
    get_ui_manager().show_warning(message)

def get_user_input(prompt: str = "请输入", mode: str = "normal") -> str:
    """全局用户输入函数"""
    return get_ui_manager().get_user_input(prompt, mode)

def show_ai_response(response: str, mode: str = "normal"):
    """全局AI响应显示函数"""
    get_ui_manager().show_ai_response(response, mode)