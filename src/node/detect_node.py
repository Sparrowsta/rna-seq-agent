from typing import Dict, Any, List
import asyncio
from ..state import AgentState
from ..tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files,
    check_fastp_availability,
    check_star_availability,
    check_hisat2_availability,
    check_featurecounts_availability
)


async def _execute_task_group(group_tasks: List[str], group_description: str) -> Dict[str, Any]:
    """执行单个任务组中的所有任务（组内串行执行）
    
    采用硬编码任务映射策略，确保执行稳定性和可预测性。
    每个检测任务都有明确的实现，避免复杂的动态调度。
    """
    print(f"🔄 开始执行{group_description}: {group_tasks}")
    
    # 硬编码任务映射 - 每个任务对应明确的检测函数
    task_mapping = {
        "analyze_fastq_data": lambda: scan_fastq_files(mode="detect", depth="detailed"),
        "assess_system_readiness": lambda: scan_system_resources(mode="detect"),
        "verify_genome_setup": lambda: scan_genome_files(mode="detect"),
        "check_fastp_availability": check_fastp_availability,
        "check_star_availability": check_star_availability,
        "check_hisat2_availability": check_hisat2_availability,
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


async def detect_node(state: AgentState) -> Dict[str, Any]:
    """Detect节点 - 直接并行执行任务组检测（硬编码策略，确保稳定性）"""
    
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