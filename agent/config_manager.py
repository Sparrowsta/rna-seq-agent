"""
配置管理模块
遵循单一职责原则：专门处理nextflow配置的读取、修改和验证
"""

import os
import json
import logging
from typing import Dict, Any, List, Optional, Union
from pathlib import Path
from dataclasses import dataclass, asdict
from copy import deepcopy

# 配置日志
logger = logging.getLogger(__name__)

# ============================================================================
# 配置数据模型 - 遵循数据类模式
# ============================================================================

@dataclass
class NextflowConfig:
    """
    Nextflow配置数据模型
    
    遵循类型安全原则：明确定义所有配置字段的类型
    """
    # 输入数据配置
    srr_ids: str = ""
    local_genome_path: str = ""
    local_gtf_path: str = ""
    download_genome_url: str = ""
    download_gtf_url: str = ""
    local_fastq_files: str = ""
    
    # 输出配置
    data: str = "./data"
    
    # 流程控制配置
    run_download_srr: bool = False
    run_download_genome: bool = False
    run_build_star_index: bool = False
    run_fastp: bool = False
    run_star_align: bool = False
    run_featurecounts: bool = False
    
    def to_dict(self) -> Dict[str, Any]:
        """转换为字典格式"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'NextflowConfig':
        """从字典创建配置对象"""
        # 过滤掉不存在的字段
        valid_fields = {field.name for field in cls.__dataclass_fields__.values()}
        filtered_data = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered_data)
    
    def validate(self) -> Dict[str, Any]:
        """
        验证配置完整性
        
        应用验证模式：集中的配置验证逻辑
        """
        errors = []
        warnings = []
        
        # 检查数据源配置
        has_fastq = bool(self.local_fastq_files.strip())
        has_srr = bool(self.srr_ids.strip())
        
        if not (has_fastq or has_srr):
            errors.append("必须配置FASTQ文件路径或SRR ID")
        
        # 检查基因组配置
        has_local_genome = bool(self.local_genome_path.strip())
        has_download_genome = bool(self.download_genome_url.strip())
        
        if not (has_local_genome or has_download_genome):
            errors.append("必须配置本地基因组路径或下载URL")
        
        # 检查GTF配置
        if has_local_genome and not self.local_gtf_path.strip():
            warnings.append("使用本地基因组时建议配置GTF文件路径")
        
        # 检查流程配置
        enabled_processes = [
            self.run_download_srr, self.run_download_genome, self.run_build_star_index,
            self.run_fastp, self.run_star_align, self.run_featurecounts
        ]
        
        if not any(enabled_processes):
            warnings.append("没有启用任何分析流程")
        
        # 检查逻辑依赖
        if self.run_star_align and not (self.run_build_star_index or has_local_genome):
            errors.append("运行STAR比对需要先构建索引或提供本地基因组")
        
        if self.run_featurecounts and not self.run_star_align:
            warnings.append("运行featureCounts通常需要先进行序列比对")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings
        }

# ============================================================================
# 配置管理器 - 遵循管理器模式
# ============================================================================

class ConfigManager:
    """
    配置管理器
    
    遵循单一职责原则：专门管理nextflow配置的CRUD操作
    """
    
    def __init__(self, config_file: str = "config/nextflow_params.json"):
        self.config_file = Path(config_file)
        self.current_config: Optional[NextflowConfig] = None
        self._ensure_config_directory()
    
    def _ensure_config_directory(self):
        """确保配置目录存在"""
        self.config_file.parent.mkdir(parents=True, exist_ok=True)
    
    def load_config(self) -> NextflowConfig:
        """
        加载配置文件
        
        应用KISS原则：简单的配置加载逻辑
        """
        try:
            if self.config_file.exists():
                with open(self.config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.current_config = NextflowConfig.from_dict(data)
                logger.info(f"配置已从 {self.config_file} 加载")
            else:
                self.current_config = NextflowConfig()
                logger.info("使用默认配置")
            
            return self.current_config
        
        except json.JSONDecodeError as e:
            logger.error(f"配置文件格式错误: {e}")
            self.current_config = NextflowConfig()
            return self.current_config
        
        except Exception as e:
            logger.error(f"加载配置文件失败: {e}")
            self.current_config = NextflowConfig()
            return self.current_config
    
    def save_config(self, config: Optional[NextflowConfig] = None) -> bool:
        """
        保存配置文件
        
        应用KISS原则：简单的配置保存逻辑
        """
        try:
            config_to_save = config or self.current_config
            if config_to_save is None:
                logger.error("没有配置可保存")
                return False
            
            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(config_to_save.to_dict(), f, indent=2, ensure_ascii=False)
            
            self.current_config = config_to_save
            logger.info(f"配置已保存到 {self.config_file}")
            return True
        
        except Exception as e:
            logger.error(f"保存配置文件失败: {e}")
            return False
    
    def update_config(self, updates: Dict[str, Any]) -> bool:
        """
        更新配置
        
        遵循DRY原则：统一的配置更新逻辑
        """
        try:
            if self.current_config is None:
                self.load_config()
            
            # 创建配置副本
            config_dict = self.current_config.to_dict()
            
            # 应用更新
            for key, value in updates.items():
                if key in config_dict:
                    config_dict[key] = value
                    logger.debug(f"更新配置: {key} = {value}")
                else:
                    logger.warning(f"忽略未知配置项: {key}")
            
            # 创建新的配置对象
            updated_config = NextflowConfig.from_dict(config_dict)
            
            # 验证配置
            validation_result = updated_config.validate()
            if not validation_result["valid"]:
                logger.error(f"配置验证失败: {validation_result['errors']}")
                return False
            
            # 保存配置
            return self.save_config(updated_config)
        
        except Exception as e:
            logger.error(f"更新配置失败: {e}")
            return False
    
    def get_config(self) -> NextflowConfig:
        """
        获取当前配置
        
        应用懒加载模式：需要时才加载配置
        """
        if self.current_config is None:
            self.load_config()
        return self.current_config
    
    def reset_config(self) -> bool:
        """
        重置为默认配置
        
        应用KISS原则：简单的重置逻辑
        """
        try:
            self.current_config = NextflowConfig()
            return self.save_config()
        except Exception as e:
            logger.error(f"重置配置失败: {e}")
            return False
    
    def export_config(self, export_path: str) -> bool:
        """
        导出配置到指定路径
        
        应用KISS原则：简单的配置导出
        """
        try:
            if self.current_config is None:
                self.load_config()
            
            export_file = Path(export_path)
            export_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(export_file, 'w', encoding='utf-8') as f:
                json.dump(self.current_config.to_dict(), f, indent=2, ensure_ascii=False)
            
            logger.info(f"配置已导出到 {export_path}")
            return True
        
        except Exception as e:
            logger.error(f"导出配置失败: {e}")
            return False
    
    def import_config(self, import_path: str) -> bool:
        """
        从指定路径导入配置
        
        应用KISS原则：简单的配置导入
        """
        try:
            import_file = Path(import_path)
            if not import_file.exists():
                logger.error(f"导入文件不存在: {import_path}")
                return False
            
            with open(import_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            imported_config = NextflowConfig.from_dict(data)
            
            # 验证导入的配置
            validation_result = imported_config.validate()
            if validation_result["warnings"]:
                logger.warning(f"导入配置存在警告: {validation_result['warnings']}")
            
            if not validation_result["valid"]:
                logger.error(f"导入配置验证失败: {validation_result['errors']}")
                return False
            
            return self.save_config(imported_config)
        
        except Exception as e:
            logger.error(f"导入配置失败: {e}")
            return False

# ============================================================================
# 配置模板管理 - 遵循模板方法模式
# ============================================================================

class ConfigTemplate:
    """
    配置模板管理器
    
    遵循模板方法模式：提供预定义的配置模板
    """
    
    @staticmethod
    def get_local_analysis_template() -> NextflowConfig:
        """本地文件分析模板"""
        return NextflowConfig(
            local_fastq_files="data/fastq/*.fastq.gz",
            local_genome_path="data/genomes/hg38/hg38.fa",
            local_gtf_path="data/genomes/hg38/hg38.gtf",
            data="./data",
            run_fastp=True,
            run_star_align=True,
            run_build_star_index=True,
            run_featurecounts=True
        )
    
    @staticmethod
    def get_srr_download_template() -> NextflowConfig:
        """SRR数据下载分析模板"""
        return NextflowConfig(
            srr_ids="SRR123456,SRR123457",
            download_genome_url="https://example.com/genome.fa.gz",
            download_gtf_url="https://example.com/annotation.gtf.gz",
            data="./data",
            run_download_srr=True,
            run_download_genome=True,
            run_build_star_index=True,
            run_fastp=True,
            run_star_align=True,
            run_featurecounts=True
        )
    
    @staticmethod
    def get_minimal_template() -> NextflowConfig:
        """最小分析模板"""
        return NextflowConfig(
            local_fastq_files="data/fastq/*.fastq.gz",
            local_genome_path="data/genomes/hg38/hg38.fa",
            data="./data",
            run_star_align=True,
            run_featurecounts=True
        )
    
    @staticmethod
    def get_comprehensive_template() -> NextflowConfig:
        """全面分析模板"""
        return NextflowConfig(
            local_fastq_files="data/fastq/*.fastq.gz",
            local_genome_path="data/genomes/hg38/hg38.fa",
            local_gtf_path="data/genomes/hg38/hg38.gtf",
            data="./data",
            run_fastp=True,
            run_build_star_index=True,
            run_star_align=True,
            run_featurecounts=True
        )
    
    @classmethod
    def list_templates(cls) -> Dict[str, str]:
        """列出所有可用模板"""
        return {
            "local": "本地文件分析模板",
            "srr": "SRR数据下载分析模板",
            "minimal": "最小分析模板",
            "comprehensive": "全面分析模板"
        }
    
    @classmethod
    def get_template(cls, template_name: str) -> Optional[NextflowConfig]:
        """根据名称获取模板"""
        templates = {
            "local": cls.get_local_analysis_template,
            "srr": cls.get_srr_download_template,
            "minimal": cls.get_minimal_template,
            "comprehensive": cls.get_comprehensive_template
        }
        
        template_func = templates.get(template_name)
        return template_func() if template_func else None

# ============================================================================
# 配置验证器 - 遵循验证器模式
# ============================================================================

class ConfigValidator:
    """
    配置验证器
    
    遵循单一职责原则：专门处理配置验证逻辑
    """
    
    @staticmethod
    def validate_file_paths(config: NextflowConfig) -> Dict[str, Any]:
        """验证文件路径"""
        results = {"valid": True, "errors": [], "warnings": []}
        
        # 检查本地文件路径
        if config.local_genome_path:
            if not os.path.exists(config.local_genome_path):
                results["errors"].append(f"基因组文件不存在: {config.local_genome_path}")
                results["valid"] = False
        
        if config.local_gtf_path:
            if not os.path.exists(config.local_gtf_path):
                results["errors"].append(f"GTF文件不存在: {config.local_gtf_path}")
                results["valid"] = False
        
        # 检查FASTQ文件路径（支持通配符）
        if config.local_fastq_files:
            import glob
            fastq_files = glob.glob(config.local_fastq_files)
            if not fastq_files:
                results["warnings"].append(f"未找到匹配的FASTQ文件: {config.local_fastq_files}")
        
        return results
    
    @staticmethod
    def validate_urls(config: NextflowConfig) -> Dict[str, Any]:
        """验证URL格式"""
        results = {"valid": True, "errors": [], "warnings": []}
        
        urls_to_check = [
            ("基因组下载URL", config.download_genome_url),
            ("GTF下载URL", config.download_gtf_url)
        ]
        
        for name, url in urls_to_check:
            if url and not (url.startswith("http://") or url.startswith("https://")):
                results["errors"].append(f"{name} 格式不正确: {url}")
                results["valid"] = False
        
        return results
    
    @staticmethod
    def validate_logic_dependencies(config: NextflowConfig) -> Dict[str, Any]:
        """验证逻辑依赖关系"""
        results = {"valid": True, "errors": [], "warnings": []}
        
        # STAR比对依赖检查
        if config.run_star_align:
            if not (config.run_build_star_index or config.local_genome_path):
                results["errors"].append("运行STAR比对需要构建索引或提供本地基因组文件")
                results["valid"] = False
        
        # featureCounts依赖检查
        if config.run_featurecounts:
            if not config.run_star_align:
                results["warnings"].append("运行featureCounts通常需要先进行STAR比对")
            
            if not (config.local_gtf_path or config.download_gtf_url):
                results["errors"].append("运行featureCounts需要GTF注释文件")
                results["valid"] = False
        
        # 数据源检查
        has_data_source = bool(config.local_fastq_files or config.srr_ids)
        if not has_data_source:
            results["errors"].append("必须指定数据源（本地FASTQ文件或SRR ID）")
            results["valid"] = False
        
        return results
    
    @classmethod
    def comprehensive_validate(cls, config: NextflowConfig) -> Dict[str, Any]:
        """综合验证"""
        all_results = {
            "valid": True,
            "errors": [],
            "warnings": []
        }
        
        # 执行各项验证
        validations = [
            cls.validate_file_paths(config),
            cls.validate_urls(config),
            cls.validate_logic_dependencies(config)
        ]
        
        # 合并结果
        for result in validations:
            if not result["valid"]:
                all_results["valid"] = False
            all_results["errors"].extend(result["errors"])
            all_results["warnings"].extend(result["warnings"])
        
        return all_results

# ============================================================================
# 全局配置管理器实例 - 应用单例模式
# ============================================================================

_global_config_manager: Optional[ConfigManager] = None

def get_config_manager() -> ConfigManager:
    """
    获取全局配置管理器实例
    
    应用单例模式：确保全局只有一个配置管理器实例
    """
    global _global_config_manager
    if _global_config_manager is None:
        _global_config_manager = ConfigManager()
    return _global_config_manager

# ============================================================================
# 便捷函数 - 遵循DRY原则
# ============================================================================

def load_current_config() -> NextflowConfig:
    """加载当前配置"""
    return get_config_manager().get_config()

def save_current_config(config: NextflowConfig) -> bool:
    """保存当前配置"""
    return get_config_manager().save_config(config)

def update_current_config(updates: Dict[str, Any]) -> bool:
    """更新当前配置"""
    return get_config_manager().update_config(updates)

def validate_current_config() -> Dict[str, Any]:
    """验证当前配置"""
    config = load_current_config()
    return ConfigValidator.comprehensive_validate(config)

def apply_template(template_name: str) -> bool:
    """应用配置模板"""
    template = ConfigTemplate.get_template(template_name)
    if template is None:
        logger.error(f"未知的模板名称: {template_name}")
        return False
    
    return save_current_config(template)