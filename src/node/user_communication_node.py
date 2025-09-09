from typing import Dict, Any
from ..state import AgentState

async def user_communication_node(state: AgentState) -> Dict[str, Any]:
    """User Communication节点 - 用户交互入口"""
    print(f"\n{'='*60}")
    print(f"🔬 RNA-seq智能分析助手 - 本地FASTQ数据分析工具")
    print(f"{'='*60}")
    print(f"")
    print(f"📋 **快速开始:**")
    print(f"   /plan                   - 🚀 开始RNA-seq分析流程")
    print(f"   /plan 使用hg19基因组    - 🎯 指定分析需求开始")
    print(f"")
    print(f"📊 **项目管理:**")
    print(f"   项目概览                - 📈 查看项目整体状态")
    print(f"   FASTQ文件查询           - 📂 浏览可用的测序数据")
    print(f"   基因组信息查询          - 🧬 检查基因组配置状态")
    print(f"   历史分析                - 📚 查看已完成的分析")
    print(f"")
    print(f"🧬 **基因组配置:**")
    print(f"   添加基因组配置          - ➕ 添加新的基因组配置")
    print(f"   添加 mm10 fasta:[url] gtf:[url]     - 🔗将额外的基因组加入基因组配置中")
    print(f"")
    print(f"⚙️ **系统命令:**")
    print(f"   /help                   - ❓ 获取详细帮助信息")
    print(f"")
    print(f"   /exit                   - 🚪 退出程序")
    print(f"")
    print(f"💡 **使用提示:**")
    print(f"   • 支持中文自然语言交互，直接描述您的分析需求")
    print(f"   • 支持逐步配置，在不使用/plan 指令下，直接输入\'使用hg19基因组\'，可单纯配置而不进入自动配置流程")
    print(f"   • 系统会自动检测FASTQ文件并智能配对")
    print(f"   • 基于Docker容器化，确保分析环境一致性")
    print(f"   • 生成标准化的Nextflow流水线，可重复执行")
    print(f"{'='*60}")
    
    # 检查并显示来自normal节点的结果
    if hasattr(state, 'query_response') and state.query_response:
        print()
        print(f"🎯 {state.query_response}")
        print()
    
    # 获取用户输入
    try:
        user_input = input("请输入: ").strip()
        
        # 定义plan等价命令
        plan_prefixes = ['/plan', '/开始分析']
        
        user_input_lower = user_input.lower()
        is_plan_command = (user_input_lower in plan_prefixes or 
                          any(user_input_lower.startswith(f"{prefix} ") for prefix in plan_prefixes))
        if user_input_lower in ['/exit', '/退出']:
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "end",
                "response": "再见！",
                "status": "normal"
            }
        elif is_plan_command:
            # 优雅的参数提取 - 处理所有plan等价命令
            plan_content = ""
            for prefix in plan_prefixes:
                if user_input_lower.startswith(prefix.lower()):
                    plan_content = user_input[len(prefix):].strip()
                    break
            
            if plan_content:
                plan_user_requirements = {"raw_input": plan_content}
                response_msg = f"进入分析计划流程...\n📝 分析需求: {plan_content}"
            else:
                plan_user_requirements = {}
                response_msg = "进入分析计划流程..."
            
            return {
                "messages": [{"role": "user", "content": user_input}], 
                "routing_decision": "plan",
                "response": response_msg,
                "user_requirements": plan_user_requirements,
                "status": "plan"
            }
        else:
            # 其他输入交给normal节点处理
            return {
                "messages": [{"role": "user", "content": user_input}],
                "routing_decision": "normal", 
                "response": "正在分析您的需求...",
                "status": "normal"
            }
            
    except KeyboardInterrupt:
        return {
            "messages": [{"role": "user", "content": "KeyboardInterrupt"}],
            "routing_decision": "end",
            "response": "用户中断退出",
            "status": "interrupted"
        }