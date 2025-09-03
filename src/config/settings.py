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
        return Path("/config").exists() and Path("/data").exists()
    
    @property
    def config_dir(self) -> Path:
        """配置文件目录"""
        if self.is_container_environment:
            return Path("/config")
        else:
            return self.project_root / "config"
    
    @property 
    def data_dir(self) -> Path:
        """数据文件目录"""
        if self.is_container_environment:
            return Path("/data")
        else:
            return self.project_root / "data"
    
    # === API配置 ===
    deepseek_api_key: str = Field(default="", description="DeepSeek API密钥")
    llm_temperature: float = Field(default=0.1, description="LLM温度参数")
    
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
        return self.config_dir / "genomes.json"
    
    @property
    def runtime_config_path(self) -> Path:
        """运行时配置文件路径"""
        return self.config_dir / "runtime_config.json"
    
    @property
    def nextflow_config_path(self) -> Path:
        """Nextflow配置文件路径"""
        return self.config_dir / "nextflow.config"
    
    @property
    def templates_dir(self) -> Path:
        """模板文件目录 - Docker环境兼容"""
        if self.is_container_environment:
            return Path("/src/templates")  # Docker中的绝对路径
        else:
            return self.project_root / "src" / "templates"  # 本地开发路径
    
    def validate_environment(self) -> tuple[bool, list[str]]:
        """验证环境配置"""
        errors = []
        
        # 检查API密钥
        if not self.deepseek_api_key:
            errors.append("DEEPSEEK_API_KEY环境变量未设置")
        
        # 检查关键目录
        if not self.config_dir.exists():
            errors.append(f"配置目录不存在: {self.config_dir}")
        
        if not self.data_dir.exists():
            errors.append(f"数据目录不存在: {self.data_dir}")
        
        return len(errors) == 0, errors