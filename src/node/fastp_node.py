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
        qc_tool = nextflow_config.get("qc_tool", "")
        align_tool = nextflow_config.get("align_tool", "")
        quant_tool = nextflow_config.get("quant_tool", "")
        
        # 将验证的工作目录路径添加到nextflow_config中
        nextflow_config["validated_work_dir"] = env_validation["work_dir"]
        
        if not sample_groups:
            return {
                "status": "error",
                "response": "❌ 错误：未找到样本信息，无法进行分析",
                "env_validation": env_validation
            }
        
        # 不在节点侧计算线程等资源；直接使用 nextflow_config 中的 fastp 相关设置
        print(f"🧹 进入FastP质控阶段（MVP：不执行下游比对/定量）...")
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
    
    # 获取历史优化参数和版本信息（实现参数迭代进化）
    current_params = getattr(state, 'fastp_current_params', {}) or {}
    current_version = getattr(state, 'fastp_version', 1)
    
    # 使用新的批次处理方法，传递历史优化参数和版本号
    print("📦 正在使用重构后的批次处理方法…")
    # 注入已有历史，便于LLM生成有上下文的建议
    history_to_pass = getattr(state, 'fastp_version_history', []) or []
    if history_to_pass:
        print(f"🧠 注入历史上下文: {len(history_to_pass)} 条记录")
    batch_result = agent.run_batch(
        sample_groups,
        nextflow_config,
        current_params,
        current_version,
        version_history=history_to_pass
    )
    
    if not batch_result.get("samples"):
        return {
            "status": "error",
            "response": f"❗ FastP批次处理失败: {batch_result.get('error', '未知错误')}",
            "results": []
        }
    
    # 提取结果和建议
    summary = batch_result.get("summary", "FastP批次处理完成")
    current_params = batch_result.get("current_params", {})
    optimized_params = batch_result.get("optimized_params", {})
    # 仅保留与当前参数不同的优化项，避免显示过时/相同值
    if optimized_params and current_params:
        optimized_params = {k: v for k, v in optimized_params.items() if current_params.get(k) != v}
    next_params = batch_result.get("next_params", {})
    reasoning = batch_result.get("reasoning", "")
    success_count = batch_result.get("success_count", 0)
    total_count = batch_result.get("total", len(sample_groups))
    version_files = batch_result.get("version_files", {})
    
    # 版本管理：创建历史记录条目和更新历史列表
    next_version = current_version + (1 if optimized_params else 0)
    version_record = {
        "version": current_version,
        "timestamp": __import__('datetime').datetime.now().isoformat(),
        "params": current_params,
        "optimized_params": optimized_params,
        "reasoning": reasoning,
        "execution_result": {
            "success_count": success_count,
            "total_count": total_count,
            "success_rate": success_count / total_count if total_count > 0 else 0
        },
        "version_files": version_files
    }
    
    # 获取当前历史记录并追加新记录
    current_history = getattr(state, 'fastp_version_history', []) or []
    updated_history = current_history + [version_record]
    
    # 打印版本信息
    print(f"✅ 批次处理完成: {success_count}/{total_count} 样本成功 [v{current_version}]")
    if version_files.get("versioned"):
        print(f"📋 参数文件已保存: {version_files['versioned']}")
    
    if reasoning:
        print(f"💡 **参数优化理由:**")
        print(f"   {reasoning}")
    
    if optimized_params:
        print(f"📊 **参数优化建议 (v{current_version} -> v{next_version}):**")
        for key, value in optimized_params.items():
            print(f"   - {key}: {value}")
    else:
        print(f"📊 **参数状态**: 稳定，无新优化建议")

    # 不再展示稳定门拦截信息
    
    return {
        "status": "qc_completed",
        "response": summary,
        # 更新state中的fastp参数（实现迭代进化）
        "fastp_prev_params": current_params,           # 本次执行前的参数（用于展示 old -> new）
        "fastp_current_params": next_params,           # 本次执行后（下次使用）的参数
        "fastp_optimized_params": optimized_params,    # 本次的优化建议（仅差异）
        "fastp_applied_updates": optimized_params,     # 本次实际应用的差异
        "fastp_next_params": {},  # 清空next_params，仅显示一次
        "fastp_version": next_version,  # 更新版本号
        "fastp_version_history": updated_history,  # 更新完整的历史记录列表
        "config_reasoning": reasoning,
        "batch_results": batch_result.get("samples", []),
        "success_rate": f"{success_count}/{total_count}"
    }


# 去除完整流水线执行，暂不实现
