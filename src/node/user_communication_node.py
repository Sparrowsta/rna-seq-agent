"""用户通信节点 - 纯数字选择入口

RNA-seq智能分析助手的主入口界面，采用纯数字选择模式。
支持执行向导、工具查询、基因组管理等功能。
"""

from typing import Dict, Any
from ..state import AgentState


async def user_communication_node(state: AgentState) -> Dict[str, Any]:
    """User Communication节点 - 用户交互入口"""
    
    # 显示入口菜单
    _display_main_menu()
    
    # 检查并显示来自normal节点的结果
    if hasattr(state, 'query_response') and state.query_response:
        print("")
        print(f"🎯 {state.query_response}")
        print("")
    
    # 获取用户输入并解析
    try:
        user_input = input("请输入选择: ").strip()
        return _parse_main_menu_input(user_input, state)
        
    except KeyboardInterrupt:
        print("\n\n👋 再见！")
        return {
            "response": "用户主动退出",
            "status": "end", 
            "routing_decision": "end"
        }
    except Exception as e:
        print(f"❌ 输入处理错误: {e}")
        return {
            "response": f"输入处理错误: {e}",
            "status": "error",
            "routing_decision": "normal"
        }


def _display_main_menu():
    """显示主入口菜单（使用新样式组件）"""
    lines = []
    lines += ["-" * 60, "🔬 RNA‑seq 智能分析助手", "-" * 60]
    lines.append("📋 请选择操作")
    lines += [
        "  1) 执行分析",
        "  2) 浏览 FASTQ 文件",
        "  3) 基因组配置（仅添加）",
        "  4) 帮助",
        "  5) 退出",
        "",
        "💡 使用提示:",
        "  • 数字选择：输入对应数字进行操作",
        "  • 自由查询：直接输入问题进行智能分析",
        "  • 支持中文自然语言交互",
        "  • 基于Docker容器化，工作目录=/data，支持相对路径",
    ]
    print("\n" + "\n".join(lines))


def _parse_main_menu_input(user_input: str, state: AgentState) -> Dict[str, Any]:
    """解析主菜单用户输入"""
    
    # 尝试解析为数字选择
    try:
        choice = int(user_input)
        return _handle_numeric_choice(choice, state)
    except ValueError:
        # 非数字输入，作为自由查询处理
        return _handle_free_query(user_input, state)


def _handle_numeric_choice(choice: int, state: AgentState) -> Dict[str, Any]:
    """处理数字选择"""
    
    if choice == 1:
        # 进入执行模式
        return _handle_execute_mode_entry(state)
    
    elif choice == 2:
        # 浏览FASTQ文件
        return {
            "response": "正在扫描FASTQ文件...",
            "input": "FASTQ文件查询",
            "status": "normal",
            "routing_decision": "normal"
        }
    
    elif choice == 3:
        # 基因组配置管理
        return _handle_genome_config_management(state)
    
    elif choice == 4:
        # 查看帮助
        return {
            "response": "正在获取帮助信息...",
            "input": "帮助",
            "status": "normal", 
            "routing_decision": "normal"
        }
    
    elif choice == 5:
        # 退出程序
        print("\n👋 感谢使用RNA-seq智能分析助手！")
        return {
            "response": "用户选择退出程序",
            "status": "end",
            "routing_decision": "end"
        }
    
    else:
        # 无效选择
        print(f"❌ 无效选择：{choice}。请输入 1-5 选择操作")
        return {
            "response": f"无效选择：{choice}",
            "status": "error",
            "routing_decision": "user_communication"  # 回到入口菜单
        }


def _handle_execute_mode_entry(state: AgentState) -> Dict[str, Any]:
    """处理执行模式入口"""
    
    print("\n" + "\n".join(["-" * 50, "🚀 执行向导", "-" * 50]))
    print("请选择执行方式:")
    print("    1) 直接开始执行 (不预填需求)")
    print("    2) 先输入执行需求 (如物种、基因组、工具偏好等)")
    print("    0) 返回上级菜单")
    
    try:
        exec_choice = input("请选择: ").strip()
        exec_choice_num = int(exec_choice)
        
        if exec_choice_num == 0:
            # 返回上级菜单
            return {
                "response": "返回主菜单",
                "status": "normal",
                "routing_decision": "user_communication"
            }
        
        elif exec_choice_num == 1:
            # 直接执行
            print("✅ 直接开始执行，进入系统检测...")
            return {
                "response": "直接开始执行RNA-seq分析",
                "status": "execute_direct",
                "routing_decision": "execute"
            }
        
        elif exec_choice_num == 2:
            # 输入执行需求
            return _handle_requirements_input(state)
        
        else:
            print(f"❌ 无效选择：{exec_choice_num}。请输入 0-2")
            return _handle_execute_mode_entry(state)  # 递归重试
            
    except ValueError:
        print(f"❌ 无效输入：{exec_choice}。请输入数字")
        return _handle_execute_mode_entry(state)  # 递归重试
    except KeyboardInterrupt:
        print("\n返回主菜单")
        return {
            "response": "用户取消执行向导",
            "status": "normal",
            "routing_decision": "user_communication"
        }


def _handle_requirements_input(state: AgentState) -> Dict[str, Any]:
    """处理执行需求输入"""
    
    print("\n💡 执行需求输入指南:")
    print("   可以描述：物种、基因组版本、测序类型、工具偏好等")
    print("   示例：使用hg38、双端测序、STAR比对、FeatureCounts定量")
    print("")
    
    try:
        requirements_text = input("请输入执行需求: ").strip()
        
        if not requirements_text:
            print("❌ 需求为空，返回执行向导")
            return _handle_execute_mode_entry(state)
        
        print(f"✅ 执行需求已记录: {requirements_text}")
        print("进入系统检测...")
        
        return {
            "response": f"已记录执行需求并开始分析: {requirements_text}",
            "input": requirements_text,
            "status": "execute_with_requirements",
            "routing_decision": "execute",
            "user_requirements": {
                "raw_input": requirements_text,
                "source": "execute_wizard"
            }
        }
        
    except KeyboardInterrupt:
        print("\n返回执行向导")
        return _handle_execute_mode_entry(state)


def _handle_genome_config_management(state: AgentState) -> Dict[str, Any]:
    """处理基因组配置管理"""
    
    # 使用state参数避免未使用警告
    _ = state
    
    # 这里先返回一个查询基因组状态的请求
    # 实际的子菜单会在normal节点中由工具调用来处理
    return {
        "response": "正在检查基因组配置状态...",
        "input": "基因组配置管理",
        "status": "normal",
        "routing_decision": "normal",
        "request_type": "genome_management"
    }


def _handle_free_query(user_input: str, state: AgentState) -> Dict[str, Any]:
    """处理自由查询输入"""
    
    # 使用state参数避免未使用警告
    _ = state
    
    print(f"🧠 正在分析您的需求: {user_input}")
    
    return {
        "response": f"正在分析您的需求: {user_input}",
        "input": user_input,
        "status": "normal",
        "routing_decision": "normal",
        # 统一对话输入格式：始终以消息列表形式提供给Agent
        "messages": [{"role": "user", "content": user_input}]
    }
