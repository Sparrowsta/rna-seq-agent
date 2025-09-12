"""视图模型定义

定义用户确认界面的结构化数据模型，分离数据与展示逻辑。
"""

from typing import Dict, Any, List, Optional, Literal
from pydantic import BaseModel, Field


class CommandHint(BaseModel):
    """命令提示信息"""
    command: str                    # 命令文本 (如 "/continue")
    description: str                # 命令描述 (如 "继续到下一步")
    icon: str = ""                  # 图标 (如 "➡️")
    available: bool = True          # 是否可用
    index: Optional[int] = None     # 数字索引 (用于纯数字选择模式)


class ParamItem(BaseModel):
    """参数项信息"""
    key: str                        # 参数键名
    value: Any                      # 参数值
    old_value: Optional[Any] = None # 旧值（用于显示变化）
    applied_optimization: bool = False  # 是否为应用的优化
    
    
class Section(BaseModel):
    """配置展示区域"""
    title: str                      # 区域标题
    icon: str = ""                  # 区域图标
    items: List[ParamItem] = Field(default_factory=list)  # 参数项列表
    visible: bool = True            # 是否显示
    
    # 三层结构
    effective: List[ParamItem] = Field(default_factory=list)      # 当前生效
    user_mods: List[ParamItem] = Field(default_factory=list)      # 用户修改
    optimizations: List[ParamItem] = Field(default_factory=list)  # 优化建议
    
    reasoning_text: Optional[str] = None  # 优化理由文本


class ResourceItem(BaseModel):
    """资源配置项"""
    process_name: str               # 进程名称
    display_name: str               # 显示名称
    cpus: str                       # CPU数量
    memory: str                     # 内存大小
    reasoning: Optional[str] = None # 配置理由


class SummaryItem(BaseModel):
    """配置摘要项"""
    key: str                        # 配置键
    value: Any                      # 配置值
    label: str                      # 显示标签
    icon: str = ""                  # 图标
    visible: bool = True            # 是否显示


class ConfirmView(BaseModel):
    """用户确认界面视图模型"""
    # 配置摘要
    summary: List[SummaryItem] = Field(default_factory=list)
    
    # 资源配置
    resources: List[ResourceItem] = Field(default_factory=list)
    
    # 参数对比区域 (FastP/STAR/FeatureCounts)
    sections: List[Section] = Field(default_factory=list)
    
    # 配置理由
    config_reasoning: Optional[str] = None
    
    # 动态命令提示
    commands: List[CommandHint] = Field(default_factory=list)
    
    # 执行进度信息
    completed_steps: List[str] = Field(default_factory=list)
    current_step: Optional[str] = None
    
    # 批次优化状态
    batch_optimization_complete: bool = False
    batch_optimizations_count: int = 0


class ConfirmDecision(BaseModel):
    """用户确认决策结果"""
    decision: str  # 改为普通字符串，不限制具体值
    payload: Dict[str, Any] = Field(default_factory=dict)
    execution_mode: Optional[str] = None  # single/optimized/batch_optimize
    modify_content: Optional[str] = None  # modify命令的具体内容
    
    # 添加决策类型判断的辅助属性
    @property
    def is_execute(self) -> bool:
        return self.decision == 'execute'
    
    @property 
    def is_continue(self) -> bool:
        return self.decision.startswith('continue')
    
    @property
    def is_workflow_step(self) -> bool:
        return self.decision in {'fastp', 'star', 'featurecounts', 'analysis'}
    
    @property
    def is_control_command(self) -> bool:
        return self.decision in {'modify', 'cancel', 'quit', 'restart'}

    # 添加需要二次输入检测的属性
    @property
    def needs_modify_content(self) -> bool:
        """检查是否需要采集修改内容"""
        return (self.decision == 'modify' and 
                self.payload.get('needs_modify_content', False))