"""
应用配置管理
集中管理所有配置项，支持环境变量和默认值
"""

import os
from pathlib import Path
from pydantic import BaseModel, Field

class Settings(BaseModel):
    """应用程序配置设置"""
    
    # === 项目路径配置 ===
    project_root: Path = Field(default_factory=lambda: Path.cwd())
    
    # === 环境检测 ===
    @property
    def is_container_environment(self) -> bool:
        """检测是否在容器环境中运行"""
        # 不使用环境变量；仅基于常见路径特征判断
        return Path("/data").exists() or Path("/.dockerenv").exists() or Path("/src").exists()
    
    
    @property 
    def data_dir(self) -> Path:
        """数据文件目录"""
        if self.is_container_environment:
            return Path("/data")
        else:
            return self.project_root / "data"
    
    # === API配置 ===
    deepseek_api_key: str = Field(default="", description="DeepSeek API密钥")
    llm_temperature: float = Field(default=0.0, description="LLM温度参数（默认0，确保确定性）")
    
    # === 系统配置 ===
    debug_mode: bool = Field(default=False, description="调试模式")
    log_level: str = Field(default="INFO", description="日志级别")
    
    def __init__(self, **kwargs):
        # 从环境变量加载配置
        env_values = {
            "deepseek_api_key": os.getenv("DEEPSEEK_API_KEY", ""),
            "debug_mode": os.getenv("DEBUG", "false").lower() == "true",
            "log_level": os.getenv("LOG_LEVEL", "INFO")
        }
        
        # 合并kwargs和环境变量
        merged_values = {**env_values, **kwargs}
        super().__init__(**merged_values)
    
    @property
    def genomes_config_path(self) -> Path:
        """基因组配置文件路径"""
        return Path("/src/genomes.json")
    
    # runtime_config_path 与 nextflow_config_path 已废弃
    
    # 模板目录已废弃
    
    def validate_environment(self) -> tuple[bool, list[str]]:
        """验证环境配置"""
        errors = []
        
        # 检查API密钥
        if not self.deepseek_api_key:
            errors.append("DEEPSEEK_API_KEY环境变量未设置")
        
        # 检查关键文件：genomes.json（新的统一位置）
        gpath = self.genomes_config_path
        if not gpath.exists():
            errors.append(f"基因组配置文件不存在: {gpath}")
        
        # 检查关键目录（config 目录已废弃，不再校验）
        if not self.data_dir.exists():
            errors.append(f"数据目录不存在: {self.data_dir}")
        
        return len(errors) == 0, errors
