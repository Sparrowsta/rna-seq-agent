from typing import Dict, Any, List
import asyncio
from ..state import AgentState, DetectResponse
from langgraph.prebuilt import create_react_agent
from langchain.tools import Tool
from ..tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files,
    check_fastp_availability,
    check_star_availability,
    check_featurecounts_availability
)
from ..core import get_shared_llm


async def _execute_task_group(group_tasks: List[str], group_description: str) -> Dict[str, Any]:
    """执行单个任务组中的所有任务（组内串行执行）"""
    print(f"🔄 开始执行{group_description}: {group_tasks}")
    
    # 任务名称到函数的映射
    task_mapping = {
        "analyze_fastq_data": lambda: scan_fastq_files(mode="detect", depth="detailed"),
        "assess_system_readiness": lambda: scan_system_resources(mode="detect"),
        "verify_genome_setup": lambda: scan_genome_files(mode="detect"),
        "check_fastp_availability": check_fastp_availability,
        "check_star_availability": check_star_availability,
        "check_featurecounts_availability": check_featurecounts_availability
    }
    
    group_results = {}
    group_errors = []
    
    # 组内串行执行任务
    for task_name in group_tasks:
        try:
            if task_name in task_mapping:
                print(f"  🎯 执行任务: {task_name}")
                result = task_mapping[task_name]()
                group_results[task_name] = result
                print(f"  ✅ {task_name} 执行完成")
            else:
                error_msg = f"未知任务: {task_name}"
                group_errors.append(error_msg)
                print(f"  ❌ {error_msg}")
        except Exception as e:
            error_msg = f"{task_name} 执行失败: {str(e)}"
            group_errors.append(error_msg)
            print(f"  ❌ {error_msg}")
    
    print(f"✅ {group_description}执行完成，成功{len(group_results)}个，失败{len(group_errors)}个")
    
    return {
        "group_description": group_description,
        "group_tasks": group_tasks,
        "results": group_results,
        "errors": group_errors,
        "success_count": len(group_results),
        "total_count": len(group_tasks)
    }


def create_detection_agent():
    """创建检测执行Agent（保留用于向后兼容）"""
    llm = get_shared_llm()
    
    # 直接定义检测工具列表
    tools = [
        Tool(
            name="analyze_fastq_data",
            func=lambda query="": scan_fastq_files(mode="detect", depth="detailed"),
            description="扫描和分析FASTQ文件。收集项目中所有FASTQ文件的信息，包括文件大小、样本配对关系、测序类型等。"
        ),
        Tool(
            name="assess_system_readiness", 
            func=lambda query="": scan_system_resources(mode="detect"),
            description="检测系统硬件资源。评估CPU核心数、内存容量、磁盘空间、系统负载等硬件信息。"
        ),
        Tool(
            name="verify_genome_setup",
            func=lambda query="": scan_genome_files(mode="detect"), 
            description="验证基因组文件配置。检查已配置基因组的FASTA文件、GTF文件、STAR索引文件的存在性和完整性。"
        ),
        Tool(
            name="check_fastp_availability",
            func=check_fastp_availability,
            description="检测fastp质控工具的可用性。在micromamba环境中测试fastp命令是否可执行。"
        ),
        Tool(
            name="check_star_availability", 
            func=check_star_availability,
            description="检测STAR比对工具的可用性。在micromamba环境中测试STAR命令是否可执行。"
        ),
        Tool(
            name="check_featurecounts_availability",
            func=check_featurecounts_availability,
            description="检测featureCounts定量工具的可用性。在micromamba环境中测试featureCounts工具是否可执行。"
        )
    ]
    
    system_prompt = """你是检测执行专家。你的任务是根据计划列表，智能执行检测任务并收集结果。

执行原则：
1. 按计划列表中的任务名称，依次调用对应的检测工具
2. 对于每个任务，只调用一次相应的工具
3. 如果检测失败，记录错误但继续执行其他任务
4. 收集所有检测结果，整合成统一的数据结构

返回格式示例：
{
  "query_results": {"检测结果按任务整理": "工具已优化输出格式"},
  "query_summary": "检测完成：基因组hg19可用，工具就绪，发现6个FASTQ文件"
}

可用的检测工具：
- analyze_fastq_data: 分析FASTQ文件
- assess_system_readiness: 检测系统资源
- verify_genome_setup: 验证基因组配置
- check_fastp_availability: 检测fastp工具
- check_star_availability: 检测STAR工具
- check_featurecounts_availability: 检测featureCounts工具

请按照计划列表执行检测，并返回JSON格式的结果。"""
    
    agent = create_react_agent(
        model=llm,
        tools=tools,
        prompt=system_prompt,
        response_format=DetectResponse
    )
    return agent


async def detect_node(state: AgentState) -> Dict[str, Any]:
    """Detect节点 - 并行执行任务组检测"""
    
    # 获取任务组信息（plan现在直接是List[List[str]]格式）
    task_groups = getattr(state, 'plan', [])
    group_descriptions = getattr(state, 'group_descriptions', [])
    execution_strategy = getattr(state, 'execution_strategy', 'parallel')
    
    # 检查是否有任务需要执行
    if not task_groups or not any(group for group in task_groups):
        return {
            "query_summary": "没有检测任务需要执行",
            "status": "prepare",
            "query_results": {}
        }
        
    # 确保有对应的组描述
    if len(group_descriptions) < len(task_groups):
        group_descriptions.extend([f"任务组{i+1}" for i in range(len(group_descriptions), len(task_groups))])
    
    print(f"🔍 Detect节点开始并行执行 {len(task_groups)} 个任务组")
    print(f"📋 执行策略: {execution_strategy}")
    
    try:
        if execution_strategy == "parallel" and len(task_groups) > 1:
            # 并行执行所有任务组
            print(f"🚀 启动 {len(task_groups)} 个任务组的并行执行")
            
            # 创建并发任务
            group_tasks = []
            for i, (group_tasks_list, group_desc) in enumerate(zip(task_groups, group_descriptions)):
                group_tasks.append(_execute_task_group(group_tasks_list, f"{group_desc}(组{i+1})"))
            
            # 使用asyncio.gather并行执行
            group_results = await asyncio.gather(*group_tasks, return_exceptions=True)
            
        else:
            # 顺序执行任务组
            print(f"🔄 顺序执行 {len(task_groups)} 个任务组")
            group_results = []
            for i, (group_tasks_list, group_desc) in enumerate(zip(task_groups, group_descriptions)):
                result = await _execute_task_group(group_tasks_list, f"{group_desc}(组{i+1})")
                group_results.append(result)
        
        # 聚合所有结果
        print(f"📊 聚合检测结果")
        all_query_results = {}
        all_errors = []
        total_success = 0
        total_tasks = 0
        
        for result in group_results:
            if isinstance(result, Exception):
                error_msg = f"任务组执行异常: {str(result)}"
                all_errors.append(error_msg)
                print(f"❌ {error_msg}")
                continue
            
            # 确保result是字典类型
            if not isinstance(result, dict):
                continue
                
            # 合并任务结果
            if 'results' in result:
                all_query_results.update(result['results'])
                total_success += result.get('success_count', 0)
                total_tasks += result.get('total_count', 0)
                
            # 收集错误
            if 'errors' in result:
                all_errors.extend(result['errors'])
        
        # 生成执行总结
        summary_parts = []
        summary_parts.append(f"并行检测完成")
        summary_parts.append(f"成功执行 {total_success}/{total_tasks} 个任务")
        
        if all_errors:
            summary_parts.append(f"发现 {len(all_errors)} 个错误")
        
        query_summary = "，".join(summary_parts)
        
        print(f"✅ 检测节点执行完成: {query_summary}")
        
        return {
            "query_summary": query_summary,
            "status": "prepare", 
            "query_results": all_query_results,
            "execution_errors": all_errors if all_errors else None,
            "execution_stats": {
                "total_groups": len(task_groups),
                "total_tasks": total_tasks,
                "successful_tasks": total_success,
                "failed_tasks": total_tasks - total_success,
                "execution_strategy": execution_strategy
            }
        }
        
    except Exception as e:
        print(f"❌ 并行检测执行失败: {str(e)}")
        import traceback
        traceback.print_exc()
        
        return {
            "query_summary": f"并行检测失败: {str(e)}",
            "status": "error",
            "query_results": {},
            "execution_errors": [str(e)]
        }