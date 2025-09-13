"""
工具配置管理
为tools.py提供配置支持，消除硬编码路径
"""

from pathlib import Path
from .settings import Settings
from .paths import PathManager

class ToolsConfig:
    """工具配置管理器"""
    
    def __init__(self, settings: Settings = None):
        self.settings = settings if settings is not None else Settings()
        self.path_manager = PathManager(self.settings)
    
    # === 数据路径配置 ===
    @property
    def project_root(self) -> Path:
        """项目根目录"""
        return self.settings.project_root
    
    @property
    def fastq_dir(self) -> Path:
        """FASTQ文件目录"""
        return self.path_manager.resolve_data_path("fastq")
    
    @property
    def results_dir(self) -> Path:
        """结果文件目录"""
        return self.path_manager.resolve_data_path("results")
    
    @property
    def reports_dir(self) -> Path:
        """报告文件目录"""
        return self.path_manager.resolve_data_path("reports")
    
    # === 配置文件路径 ===
    @property
    def genomes_config_path(self) -> Path:
        """基因组配置文件路径"""
        return self.settings.genomes_config_path
    
    # 运行时与 nextflow 配置文件路径已废弃，不再暴露属性
    
    # === 索引目录路径 ===
    def get_star_index_dir(self, fasta_path: Path) -> Path:
        """获取STAR索引目录"""
        return fasta_path.parent / "star_index"
    
    def get_hisat2_index_dir(self, fasta_path: Path) -> Path:
        """获取HISAT2索引目录"""
        return fasta_path.parent / "hisat2_index"
    
    # === 确保目录存在 ===
    def ensure_directories(self):
        """确保所有必要目录存在"""
        # 注意：config 目录已废弃，不再确保其存在
        directories = [
            self.settings.data_dir,
            self.fastq_dir,
            self.results_dir,
            self.reports_dir
        ]
        
        for directory in directories:
            self.path_manager.ensure_directory(directory)

# 全局工具配置实例
_global_tools_config = None

def get_tools_config() -> ToolsConfig:
    """获取全局工具配置实例"""
    global _global_tools_config
    if _global_tools_config is None:
        _global_tools_config = ToolsConfig()
    return _global_tools_config
