"""参数差异对比工具

提供通用的字典扁平化和三层对比（Base/Mods/Opt）功能，
减少重复代码并提升可测试性。
"""

from typing import Dict, Any, Optional, List
from pydantic import BaseModel


class ParamDiff(BaseModel):
    """参数差异对比结果"""
    effective: Dict[str, Any]           # 当前生效值
    user_modifications: Dict[str, Any]  # 用户修改项
    optimizations: Dict[str, Any]       # 优化应用项
    
    
def flatten_dict(d: Dict[str, Any], parent: str = "", sep: str = ".") -> Dict[str, Any]:
    """
    递归扁平化字典结构
    
    Args:
        d: 待扁平化的字典
        parent: 父级键名
        sep: 分隔符
        
    Returns:
        扁平化后的字典
        
    Examples:
        >>> flatten_dict({"a": {"b": 1, "c": 2}})
        {"a.b": 1, "a.c": 2}
    """
    if not isinstance(d, dict):
        return {} if d is None else {parent: d} if parent else {}
        
    out: Dict[str, Any] = {}
    for k, v in d.items():
        key = f"{parent}{sep}{k}" if parent else str(k)
        if isinstance(v, dict):
            out.update(flatten_dict(v, key, sep))
        else:
            out[key] = v
    return out


def diff_base_mods_opt(
    base: Dict[str, Any],
    effective: Dict[str, Any],
    user_mods: Optional[Dict[str, Any]] = None,
    opt: Optional[Dict[str, Any]] = None
) -> ParamDiff:
    """
    执行Base/Mods/Opt三层参数对比
    
    Args:
        base: 基线参数（首次生成的默认值）
        effective: 当前生效参数
        user_mods: 用户修改记录（可选）
        opt: 优化建议应用记录（可选）
        
    Returns:
        ParamDiff对象，包含三层对比结果
    """
    # 确保输入都是字典
    base = base or {}
    effective = effective or {}
    user_mods = user_mods or {}
    opt = opt or {}
    
    # 扁平化处理
    flattened_effective = flatten_dict(effective)
    
    # 构建用户修改项：优先使用明确的user_mods，否则与base对比推断
    user_modifications: Dict[str, Any] = {}
    if user_mods:
        user_modifications = dict(user_mods)
    else:
        # 回退方式：通过base与effective差异推断用户修改
        flattened_base = flatten_dict(base)
        for key, eff_value in flattened_effective.items():
            base_value = flattened_base.get(key)
            if base_value != eff_value:
                user_modifications[key] = eff_value
    
    # 构建优化项
    optimizations = dict(opt) if opt else {}
    
    return ParamDiff(
        effective=flattened_effective,
        user_modifications=user_modifications,
        optimizations=optimizations
    )


def extract_user_modifications_from_history(
    modification_history: List[Dict[str, Any]], 
    tool_name: str
) -> Dict[str, Any]:
    """
    从修改历史中提取指定工具的最近用户修改
    
    Args:
        modification_history: 修改历史记录列表
        tool_name: 工具名称 (fastp/star/featurecounts)
        
    Returns:
        最近的用户修改字典
    """
    if not modification_history:
        return {}
    
    try:
        last_record = modification_history[-1] or {}
        changes = last_record.get('changes', {}) or {}
        tool_changes = changes.get(tool_name, {}) or {}
        return dict(tool_changes) if isinstance(tool_changes, dict) else {}
    except (IndexError, KeyError, AttributeError):
        return {}


def extract_optimization_from_history(
    params_history: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """
    从参数历史中提取最新的优化应用记录
    
    Args:
        params_history: 参数历史记录列表
        
    Returns:
        最新的优化应用字典
    """
    if not params_history:
        return {}
        
    try:
        last_entry = params_history[-1] or {}
        optimization_applied = last_entry.get('optimization_applied', {}) or {}
        return dict(optimization_applied) if isinstance(optimization_applied, dict) else {}
    except (IndexError, KeyError, AttributeError):
        return {}


def filter_ignored_keys(flattened_dict: Dict[str, Any], ignore_keys: set) -> Dict[str, Any]:
    """
    过滤掉指定的忽略键
    
    Args:
        flattened_dict: 扁平化字典
        ignore_keys: 要忽略的键集合
        
    Returns:
        过滤后的字典
    """
    return {
        k: v for k, v in flattened_dict.items() 
        if k.split('.')[0] not in ignore_keys and k not in ignore_keys
    }
