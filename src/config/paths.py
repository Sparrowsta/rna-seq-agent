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
    
    def resolve_config_path(self, relative_path: Union[str, Path]) -> Path:
        """解析配置文件路径
        
        Args:
            relative_path: 相对于配置目录的路径
            
        Returns:
            完整的绝对路径
        """
        return self.settings.config_dir / relative_path
    
    def resolve_fastq_path(self, filename: str) -> Path:
        """解析FASTQ文件路径"""
        return self.resolve_data_path("fastq") / filename
    
    def resolve_results_path(self, relative_path: Union[str, Path]) -> Path:
        """解析结果文件路径"""
        return self.resolve_data_path("results") / relative_path
    
    def resolve_reports_path(self, relative_path: Union[str, Path]) -> Path:
        """解析报告文件路径"""
        return self.resolve_data_path("reports") / relative_path
    
    def ensure_directory(self, path: Path) -> Path:
        """确保目录存在，如果不存在则创建
        
        Args:
            path: 目录路径
            
        Returns:
            创建后的目录路径
        """
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    def ensure_parent_directory(self, file_path: Path) -> Path:
        """确保文件的父目录存在
        
        Args:
            file_path: 文件路径
            
        Returns:
            文件路径（父目录已确保存在）
        """
        file_path.parent.mkdir(parents=True, exist_ok=True)
        return file_path
    
    def get_relative_to_project(self, absolute_path: Path) -> Path:
        """获取相对于项目根目录的相对路径
        
        Args:
            absolute_path: 绝对路径
            
        Returns:
            相对于项目根目录的路径
        """
        try:
            return absolute_path.relative_to(self.settings.project_root)
        except ValueError:
            # 如果路径不在项目根目录下，返回绝对路径
            return absolute_path
    
    def is_fastq_file(self, file_path: Path) -> bool:
        """判断是否为FASTQ文件"""
        valid_extensions = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}
        
        # 处理.gz压缩文件
        if file_path.suffix == '.gz':
            stem_with_ext = file_path.stem
            return any(stem_with_ext.endswith(ext) for ext in ['.fastq', '.fq'])
        
        return file_path.suffix in valid_extensions