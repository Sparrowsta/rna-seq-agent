"""
FastP Node - 专门的FastP质控处理节点
支持智能参数配置和结果分析，可从user_confirm_node直接路由
"""

import json
import os
from typing import Dict, Any, List
from pathlib import Path

from ..state import AgentState
from ..agents.fastp_agent import FastpAgent, FastpConfig
from ..config.settings import Settings


def validate_nextflow_environment(settings: Settings) -> Dict[str, Any]:
    """验证Nextflow工作环境和变量设置
    
    Args:
        settings: 应用配置实例
        
    Returns:
        包含验证结果的字典
    """
    validation_result = {
        "status": "success",
        "errors": [],
        "warnings": [],
        "work_dir": "",
        "env_vars": {}
    }
    
    # 获取NXF_WORK环境变量
    nxf_work = os.getenv("NXF_WORK")
    nxf_home = os.getenv("NXF_HOME")
    
    if nxf_work:
        work_dir = Path(nxf_work)
        validation_result["work_dir"] = str(work_dir)
        validation_result["env_vars"]["NXF_WORK"] = nxf_work
        
        # 检查工作目录是否存在，不存在则创建
        try:
            work_dir.mkdir(parents=True, exist_ok=True)
            print(f"📁 Nextflow工作目录验证成功: {work_dir}")
        except PermissionError:
            error_msg = f"无权限创建Nextflow工作目录: {work_dir}"
            validation_result["errors"].append(error_msg)
            validation_result["status"] = "error"
        except Exception as e:
            error_msg = f"创建Nextflow工作目录失败: {work_dir}, 错误: {str(e)}"
            validation_result["errors"].append(error_msg)
            validation_result["status"] = "error"
    else:
        # 使用默认工作目录
        default_work_dir = settings.data_dir / "work"
        validation_result["work_dir"] = str(default_work_dir)
        validation_result["warnings"].append("NXF_WORK环境变量未设置，使用默认路径")
        
        try:
            default_work_dir.mkdir(parents=True, exist_ok=True)
            print(f"📁 使用默认Nextflow工作目录: {default_work_dir}")
        except Exception as e:
            error_msg = f"创建默认工作目录失败: {default_work_dir}, 错误: {str(e)}"
            validation_result["errors"].append(error_msg)
            validation_result["status"] = "error"
    
    # 记录其他相关环境变量
    if nxf_home:
        validation_result["env_vars"]["NXF_HOME"] = nxf_home
    
    return validation_result


async def fastp_node(state: AgentState) -> Dict[str, Any]:
    """
    统一执行节点 - 处理所有分析任务
    
    从user_confirm_node路由而来，根据用户确认的配置执行相应的分析流程：
    - 质控任务：使用FastP进行质控
    - 完整流水线：执行质控+比对+定量的完整分析
    """
    print(f"🚀 开始分析任务处理...")
    
    # 验证Nextflow环境
    settings = Settings()
    env_validation = validate_nextflow_environment(settings)
    
    if env_validation["status"] == "error":
        error_msg = "Nextflow环境验证失败: " + "; ".join(env_validation["errors"])
        print(f"❌ {error_msg}")
        return {
            "status": "error",
            "response": f"❌ {error_msg}",
            "env_validation": env_validation
        }
    
    # 显示环境验证警告
    if env_validation["warnings"]:
        for warning in env_validation["warnings"]:
            print(f"⚠️ {warning}")
    
    try:
        # 获取配置信息
        nextflow_config = state.nextflow_config or {}
        sample_groups = nextflow_config.get("sample_groups", [])
        
        # 将验证的工作目录路径添加到nextflow_config中
        nextflow_config["validated_work_dir"] = env_validation["work_dir"]
        
        if not sample_groups:
            return {
                "status": "error",
                "response": "❌ 错误：未找到样本信息，无法进行分析",
                "env_validation": env_validation
            }
        
        # 不在节点侧计算线程等资源；直接使用 nextflow_config 中的 fastp 相关设置
        print(f"🧹 进入FastP质控阶段（不执行下游比对/定量）...")
        result = await _execute_qc_only(state, nextflow_config, sample_groups)
        
        # 将环境验证信息添加到结果中
        result["env_validation"] = env_validation
        return result
            
    except Exception as e:
        error_msg = f"分析任务处理出错: {str(e)}"
        print(f"❌ {error_msg}")
        import traceback
        traceback.print_exc()
        
        return {
            "status": "error",
            "response": f"❌ {error_msg}",
            "results": []
        }


async def _execute_qc_only(state: AgentState, nextflow_config: Dict[str, Any], sample_groups: List[Dict]) -> Dict[str, Any]:
    """使用重构后的FastpAgent进行批次质量控制"""
    from ..agents.fastp_agent import FastpAgent
    
    # 创建FastP Agent
    agent = FastpAgent()
    
    # 获取当前FastP参数（使用简化的单一参数集）
    current_params = getattr(state, 'fastp_params', {}) or {}
    
    # 版本号仅用于历史记录，不再作为主要版本管理机制
    current_version = len(getattr(state, 'fastp_params_history', [])) + 1
    
    # 使用新的批次处理方法
    print("📦 正在使用重构后的批次处理方法…")
    
    # 传递历史信息供参考
    params_history = getattr(state, 'fastp_params_history', []) or []
    if params_history:
        print(f"🧠 参考历史记录: {len(params_history)} 次修改")
    
    batch_result = agent.run_batch(
        sample_groups,
        nextflow_config,
        current_params,
        current_version,
        version_history=params_history  # 传递简化的历史记录
    )
    
    if not batch_result.get("samples"):
        return {
            "status": "error",
            "response": f"❗ FastP批次处理失败: {batch_result.get('error', '未知错误')}",
            "results": []
        }
    
    # 提取结果
    summary = batch_result.get("summary", "FastP批次处理完成")
    execution_params = batch_result.get("current_params", {}) or current_params
    optimized_params = batch_result.get("optimized_params", {})
    reasoning = batch_result.get("reasoning", "")
    success_count = batch_result.get("success_count", 0)
    total_count = batch_result.get("total", len(sample_groups))
    version_files = batch_result.get("version_files", {})
    
    # 生成下次运行的参数（基于优化建议）
    next_run_params = execution_params.copy()
    if optimized_params:
        # 如果有优化建议，合并到参数中
        next_run_params.update(optimized_params)
        params_updated = True
    else:
        params_updated = False
    
    # 记录执行历史
    history_record = {
        "timestamp": __import__('datetime').datetime.now().isoformat(),
        "params_used": execution_params,  # 本次实际使用的参数
        "optimization_applied": optimized_params if params_updated else {},
        "reasoning": reasoning,
        "execution_result": {
            "success_count": success_count,
            "total_count": total_count,
            "success_rate": success_count / total_count if total_count > 0 else 0
        }
    }
    
    # 更新历史记录
    current_history = getattr(state, 'fastp_params_history', []) or []
    updated_history = current_history + [history_record]
    
    # 打印执行结果
    print(f"✅ 批次处理完成: {success_count}/{total_count} 样本成功")
    if version_files.get("versioned"):
        print(f"📋 参数文件已保存: {version_files['versioned']}")
    
    if reasoning:
        print(f"💡 **优化分析:**")
        print(f"   {reasoning}")
    
    # 根据执行模式决定行为
    execution_mode = getattr(state, 'execution_mode', 'single')
    
    if execution_mode == 'optimized' and params_updated:
        # 优化模式：直接应用优化参数准备下次运行
        print(f"\n📊 **参数优化应用:**")
        for key, value in optimized_params.items():
            old_val = execution_params.get(key, "未设置")
            print(f"   - {key}: {old_val} → {value}")
        print(f"⚡ 优化参数已应用，准备下次迭代")
        
        return {
            "status": "qc_completed",
            "response": f"{summary}\n✨ 优化参数已应用，可继续执行优化",
            # 直接更新参数为下次运行准备
            "fastp_params": next_run_params,
            "fastp_params_history": updated_history,
            "fastp_optimized_suggestions": optimized_params,
            "config_reasoning": reasoning,
            "batch_results": batch_result.get("samples", []),
            "success_rate": f"{success_count}/{total_count}"
        }
    else:
        # 单次模式或无优化：保持当前参数
        if params_updated:
            print(f"\n📊 **优化建议（供参考）:**")
            for key, value in optimized_params.items():
                old_val = execution_params.get(key, "未设置")
                print(f"   - {key}: {old_val} → {value}")
            print(f"💡 提示: 使用 /execute_opt 进行优化迭代")
        else:
            print(f"📊 **参数状态**: 当前参数运行良好")
        
        return {
            "status": "qc_completed",
            "response": summary,
            # 单次模式：保持原参数不变
            "fastp_params": execution_params,
            "fastp_params_history": updated_history,
            "fastp_optimized_suggestions": optimized_params,
            "config_reasoning": reasoning,
            "batch_results": batch_result.get("samples", []),
            "success_rate": f"{success_count}/{total_count}"
        }

