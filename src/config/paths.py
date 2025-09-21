"""
路径管理工具
提供统一的路径解析和管理功能
"""

from pathlib import Path
from typing import Union
from .settings import Settings

class PathManager:
    """路径管理器 - 统一处理项目中的所有路径操作"""
    
    def __init__(self, settings: Settings):
        self.settings = settings
    
    def resolve_data_path(self, relative_path: Union[str, Path]) -> Path:
        """解析数据文件路径
        
        Args:
            relative_path: 相对于数据目录的路径
            
        Returns:
            完整的绝对路径
        """
        return self.settings.data_dir / relative_path
    
    
