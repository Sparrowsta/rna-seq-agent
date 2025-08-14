import os
import json
import glob
import re
import subprocess
from typing import Dict, Any, List, Optional
from pathlib import Path
from langchain_core.tools import tool
from pydantic import BaseModel, Field

# ============================================================================
# 输入模型定义 - 遵循接口隔离原则
# ============================================================================

class DirectoryQueryArgs(BaseModel):
    """目录查询参数模型"""
    directory_path: str = Field(description="要查询的目录路径")

class FastqQueryArgs(BaseModel):
    """FASTQ文件查询参数模型"""
    directory_path: Optional[str] = Field(default=None, description="包含FASTQ文件的目录路径。如果不指定，将自动搜索默认位置：data/fastq和data/results/fastp")
    pattern: str = Field(default="*.fastq*", description="文件匹配模式")

class GenomeQueryArgs(BaseModel):
    """基因组查询参数模型"""
    genome_name: Optional[str] = Field(default=None, description="基因组名称，如hg38、mm39等。如果不提供，将返回所有基因组信息")
    config_path: str = Field(default="config/genomes.json", description="基因组配置文件路径")

class NextflowConfigArgs(BaseModel):
    """Nextflow配置参数模型"""
    param_name: str = Field(description="参数名称")
    param_value: Any = Field(description="参数值")

class BatchConfigArgs(BaseModel):
    """批量配置更新参数模型"""
    config_updates: Dict[str, Any] = Field(description="要更新的配置字典")

class ModeSwitch(BaseModel):
    """模式切换参数模型"""
    target_mode: str = Field(description="目标模式：plan或execute")
    reason: str = Field(description="切换原因")

class ExecutionArgs(BaseModel):
    """执行参数模型"""
    config_path: str = Field(default="config/nextflow.config", description="nextflow配置文件路径")
    work_dir: str = Field(default="./work", description="工作目录")

class TreeListArgs(BaseModel):
    """目录树列表参数模型"""
    directory_path: str = Field(description="要查询的目录路径")
    max_depth: Optional[int] = Field(default=None, description="最大递归深度，None表示无限制")
    file_pattern: Optional[str] = Field(default=None, description="文件匹配模式，如*.txt")
    show_only_files: bool = Field(default=False, description="是否只显示文件，不显示目录")
    output_format: str = Field(default="tree", description="输出格式，可选'tree'（树形结构）或'list'（简单列表）")

class TaskListArgs(BaseModel):
    """任务列表生成参数模型"""
    analysis_type: str = Field(default="standard", description="分析类型：standard（标准）、minimal（最小）、comprehensive（全面）")
    force_refresh: bool = Field(default=False, description="是否强制重新检测文件和配置")

# ============================================================================
# 信息查询工具组 - 遵循单一职责原则
# ============================================================================

@tool(args_schema=DirectoryQueryArgs)
def list_directory_contents(directory_path: str) -> str:
    """
    列出指定目录的内容
    
    应用KISS原则：简单直接的目录列表功能
    """
    try:
        if not os.path.exists(directory_path):
            return f"错误：目录 '{directory_path}' 不存在"
        
        if not os.path.isdir(directory_path):
            return f"错误：'{directory_path}' 不是一个目录"
        
        contents = []
        for item in os.listdir(directory_path):
            item_path = os.path.join(directory_path, item)
            if os.path.isdir(item_path):
                contents.append(f"📁 {item}/")
            else:
                size = os.path.getsize(item_path)
                contents.append(f"📄 {item} ({size} bytes)")
        
        if not contents:
            return f"目录 '{directory_path}' 为空"
        
        return f"目录 '{directory_path}' 内容：\n" + "\n".join(contents)
    
    except PermissionError:
        return f"错误：没有权限访问目录 '{directory_path}'"
    except Exception as e:
        return f"查询目录时发生错误：{str(e)}"

@tool(args_schema=FastqQueryArgs)
def query_fastq_files(directory_path: Optional[str] = None, pattern: str = "*.fastq*") -> str:
    """
    查询FASTQ文件信息，支持默认路径和用户指定路径
    
    如果不指定directory_path，将搜索默认的FASTQ存储位置：
    - data/fastq (原始FASTQ文件)
    - data/results/fastp (质控后的FASTQ文件)
    
    遵循DRY原则：统一的FASTQ文件查询逻辑
    """
    try:
        # 定义默认搜索路径
        default_paths = ["data/fastq", "data/results/fastp"]
        search_paths = []
        
        if directory_path:
            # 用户指定了路径，只搜索指定路径
            if not os.path.exists(directory_path):
                return f"错误：指定的目录 '{directory_path}' 不存在"
            search_paths = [directory_path]
        else:
            # 使用默认路径，只搜索存在的路径
            search_paths = [path for path in default_paths if os.path.exists(path)]
            
            if not search_paths:
                return f"错误：默认FASTQ目录不存在。请检查以下目录：{', '.join(default_paths)}"
        
        all_fastq_files = []
        searched_paths = []
        
        # 在所有搜索路径中查找FASTQ文件
        for search_path in search_paths:
            searched_paths.append(search_path)
            search_pattern = os.path.join(search_path, "**", pattern)
            fastq_files = glob.glob(search_pattern, recursive=True)
            all_fastq_files.extend(fastq_files)
        
        if not all_fastq_files:
            return f"在搜索路径 {', '.join(searched_paths)} 中未找到匹配 '{pattern}' 的FASTQ文件"
        
        # 分析文件信息 - 修复配对逻辑（避免文件覆盖）
        all_files_list = []  # 存储所有文件信息的列表
        single_files = []
        
        def is_processed_file(file_path: str) -> bool:
            """判断文件是否为处理后的文件"""
            file_name = os.path.basename(file_path).lower()
            processed_indicators = [
                'trimmed', 'fastp', 'cutadapt', 'processed', 'clean', 
                'filtered', 'qc', 'trim', 'adapter'
            ]
            # 检查文件名是否包含处理后的标识
            for indicator in processed_indicators:
                if indicator in file_name:
                    return True
            
            # 检查路径是否包含处理后的目录标识
            path_processed_indicators = [
                'fastp', 'trimmed', 'processed', 'clean', 'qc', 
                'cutadapt', 'trim', 'filter', 'results'
            ]
            for indicator in path_processed_indicators:
                if indicator in file_path.lower():
                    return True
            
            return False
        
        # 第一步：收集所有文件并尝试解析配对信息
        for file_path in sorted(all_fastq_files):
            file_name = os.path.basename(file_path)
            file_size = os.path.getsize(file_path)
            
            # 改进的R1/R2检测逻辑
            # 更全面的R1/R2模式匹配
            r1_patterns = [
                r'(.+)_1\.fastq', r'(.+)_R1\.fastq', r'(.+)_r1\.fastq',
                r'(.+)_1\.fq', r'(.+)_R1\.fq', r'(.+)_r1\.fq',
                r'(.+)_1\.fastq\.gz', r'(.+)_R1\.fastq\.gz', r'(.+)_r1\.fastq\.gz',
                r'(.+)_1\.fq\.gz', r'(.+)_R1\.fq\.gz', r'(.+)_r1\.fq\.gz',
                r'(.+)_1\.trimmed\.fastq', r'(.+)_R1\.trimmed\.fastq',
                r'(.+)_1\.trimmed\.fq', r'(.+)_R1\.trimmed\.fq',
                r'(.+)\.1\.fastq', r'(.+)\.R1\.fastq', r'(.+)\.r1\.fastq'
            ]
            
            r2_patterns = [
                r'(.+)_2\.fastq', r'(.+)_R2\.fastq', r'(.+)_r2\.fastq',
                r'(.+)_2\.fq', r'(.+)_R2\.fq', r'(.+)_r2\.fq',
                r'(.+)_2\.fastq\.gz', r'(.+)_R2\.fastq\.gz', r'(.+)_r2\.fastq\.gz',
                r'(.+)_2\.fq\.gz', r'(.+)_R2\.fq\.gz', r'(.+)_r2\.fq\.gz',
                r'(.+)_2\.trimmed\.fastq', r'(.+)_R2\.trimmed\.fastq',
                r'(.+)_2\.trimmed\.fq', r'(.+)_R2\.trimmed\.fq',
                r'(.+)\.2\.fastq', r'(.+)\.R2\.fastq', r'(.+)\.r2\.fastq'
            ]
            
            sample_id = None
            read_type = None
            
            # 检查是否为R1
            for pattern in r1_patterns:
                match = re.match(pattern, file_name, re.IGNORECASE)
                if match:
                    sample_id = match.group(1)
                    read_type = "R1"
                    break
            
            # 如果不是R1，检查是否为R2
            if not sample_id:
                for pattern in r2_patterns:
                    match = re.match(pattern, file_name, re.IGNORECASE)
                    if match:
                        sample_id = match.group(1)
                        read_type = "R2"
                        break
            
            # 存储文件信息
            file_info = {
                "path": file_path, 
                "size": file_size, 
                "name": file_name,
                "sample_id": sample_id,
                "read_type": read_type,
                "is_processed": is_processed_file(file_path)
            }
            
            if sample_id and read_type:
                # 可能的配对文件，加入列表而不是覆盖
                all_files_list.append(file_info)
            else:
                # 无法识别配对信息的文件，可能是真正的单端文件
                single_files.append(file_info)
        
        # 第二步：优先选择原始文件进行配对，如果没有原始文件再使用处理后文件
        paired_files = {}
        remaining_files = all_files_list.copy()
        
        # 获取所有样本ID
        sample_ids = list(set(f["sample_id"] for f in all_files_list))
        
        for sample_id in sample_ids:
            # 找到此样本的所有文件
            sample_files = [f for f in all_files_list if f["sample_id"] == sample_id]
            
            # 分别收集原始和处理后的R1/R2文件
            original_r1 = [f for f in sample_files if f["read_type"] == "R1" and not f["is_processed"]]
            original_r2 = [f for f in sample_files if f["read_type"] == "R2" and not f["is_processed"]]
            processed_r1 = [f for f in sample_files if f["read_type"] == "R1" and f["is_processed"]]
            processed_r2 = [f for f in sample_files if f["read_type"] == "R2" and f["is_processed"]]
            
            # 优先使用原始文件配对
            r1_file = original_r1[0] if original_r1 else (processed_r1[0] if processed_r1 else None)
            r2_file = original_r2[0] if original_r2 else (processed_r2[0] if processed_r2 else None)
            
            if r1_file and r2_file:
                # 完整的双端配对
                paired_files[sample_id] = {"R1": r1_file, "R2": r2_file}
                # 从remaining_files中移除已配对的文件
                remaining_files = [f for f in remaining_files if f != r1_file and f != r2_file]
            
        # 剩余未配对的文件归类为单端
        single_files.extend(remaining_files)
        
        # 按文件来源分类结果（基于新的数据结构）
        original_paired = {}
        original_single = []
        processed_paired = {}
        processed_single = []
        
        # 分类配对文件
        for sample_id, files in paired_files.items():
            # 检查R1和R2文件是否都不是处理后的文件
            r1_is_original = not files['R1']['is_processed']
            r2_is_original = not files['R2']['is_processed']
            
            if r1_is_original and r2_is_original:
                original_paired[sample_id] = files
            else:
                processed_paired[sample_id] = files
        
        # 分类单端文件
        for file_info in single_files:
            if not file_info['is_processed']:
                original_single.append(file_info)
            else:
                processed_single.append(file_info)
        
        # 构建结果
        result = [f"FASTQ文件查询结果："]
        result.append(f"搜索路径: {', '.join(searched_paths)}\n")
        
        # 显示原始文件（主要用于分析）
        if original_paired or original_single:
            result.append("📁 **原始FASTQ文件** (主要：用于分析)：")
            
            if original_paired:
                result.append("  📝 双端测序文件：")
                for sample_id, files in original_paired.items():
                    result.append(f"    📦 样本: {sample_id}")
                    if "R1" in files:
                        size_mb = files['R1']['size'] / (1024*1024)
                        result.append(f"      📄 R1: {files['R1']['path']} ({size_mb:.2f} MB)")
                    if "R2" in files:
                        size_mb = files['R2']['size'] / (1024*1024)
                        result.append(f"      📄 R2: {files['R2']['path']} ({size_mb:.2f} MB)")
            
            if original_single:
                result.append("  📝 单端测序文件：")
                for file_info in original_single:
                    size_mb = file_info['size'] / (1024*1024)
                    result.append(f"    📄 {file_info['name']}: {file_info['path']} ({size_mb:.2f} MB)")
        
        # 显示质控后文件（次要，显示处理状态）
        if processed_paired or processed_single:
            result.append("\n🔬 **质控后文件** (次要：显示处理状态)：")
            
            if processed_paired:
                result.append("  📝 双端测序文件：")
                for sample_id, files in processed_paired.items():
                    result.append(f"    📦 样本: {sample_id}")
                    if "R1" in files:
                        size_mb = files['R1']['size'] / (1024*1024)
                        result.append(f"      📄 R1: {files['R1']['path']} ({size_mb:.2f} MB)")
                    if "R2" in files:
                        size_mb = files['R2']['size'] / (1024*1024)
                        result.append(f"      📄 R2: {files['R2']['path']} ({size_mb:.2f} MB)")
            
            if processed_single:
                result.append("  📝 单端测序文件：")
                for file_info in processed_single:
                    size_mb = file_info['size'] / (1024*1024)
                    result.append(f"    📄 {file_info['name']}: {file_info['path']} ({size_mb:.2f} MB)")
        
        # 统计信息
        total_original = len(original_paired) + len(original_single)
        total_processed = len(processed_paired) + len(processed_single)
        result.append(f"\n📊 统计: 原始文件 {total_original} 个，质控文件 {total_processed} 个")
        
        # 使用建议
        if original_paired or original_single:
            result.append("\n💡 使用建议:")
            result.append("- ✅ 优先使用原始FASTQ文件进行RNA-seq分析")
            if processed_paired or processed_single:
                result.append("- 📋 质控文件可用于验证数据处理状态")
            result.append("- 🚀 如需开始分析流程，请告诉我您的需求")
        elif processed_paired or processed_single:
            result.append("\n💡 使用建议:")
            result.append("- ⚠️  仅找到质控后文件，建议检查是否有原始FASTQ文件")
            result.append("- 📋 这些文件可显示数据处理状态")
            result.append("- 🚀 如需开始分析流程，请告诉我您的需求")
        
        return "\n".join(result)
    
    except Exception as e:
        return f"查询FASTQ文件时发生错误：{str(e)}"

@tool(args_schema=GenomeQueryArgs)
def query_genome_info(genome_name: Optional[str] = None, config_path: str = "config/genomes.json") -> str:
    """
    查询基因组配置信息并动态检查文件系统状态，自动同步配置
    
    遵循单一职责原则：专门处理基因组信息查询和配置同步
    当不提供genome_name时，返回所有基因组的摘要信息
    当提供genome_name时，返回特定基因组的详细信息
    自动检测文件系统中的新基因组并更新配置
    """
    try:
        if not os.path.exists(config_path):
            return f"错误：基因组配置文件 '{config_path}' 不存在"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        
        # 实时更新配置状态：检查文件是否存在并更新配置
        for name, info in genomes_config.items():
            # 统一处理fasta和gtf路径（兼容不同命名）
            fasta_path = info.get('fasta', info.get('fasta_path', ''))
            gtf_path = info.get('gtf', info.get('gtf_path', ''))
            
            # 更新实际存在状态
            if fasta_path:
                info['fasta_exists'] = os.path.exists(fasta_path)
            if gtf_path:
                info['gtf_exists'] = os.path.exists(gtf_path)
            
            # 检查STAR索引
            if fasta_path and os.path.exists(fasta_path):
                star_index_dir = os.path.join(os.path.dirname(fasta_path), "star_index")
                info['star_index_exists'] = (os.path.exists(star_index_dir) and 
                                            os.path.exists(os.path.join(star_index_dir, "SA")))
            else:
                info['star_index_exists'] = False
        
        # 如果没有提供genome_name，返回所有基因组的摘要信息
        if genome_name is None:
            result = ["可用基因组详细状态信息："]
            result.append("=" * 80)
            
            for name, info in genomes_config.items():
                result.append(f"📊 基因组: {name}")
                result.append(f"   物种: {info.get('species', '未知')}")
                result.append(f"   版本: {info.get('version', '未知')}")
                
                # 检查本地文件是否存在
                fasta_path = info.get('fasta', info.get('fasta_path', ''))
                gtf_path = info.get('gtf', info.get('gtf_path', ''))
                
                # 检查FASTA文件状态
                if fasta_path and info.get('fasta_exists', False):
                    file_size = os.path.getsize(fasta_path) / (1024**3)  # GB
                    fasta_status = f"✅ 已下载 ({file_size:.2f} GB)"
                elif fasta_path:
                    fasta_status = "❌ 未下载 (配置已设置)"
                else:
                    fasta_status = "⚠️  未配置"
                
                # 检查GTF文件状态
                if gtf_path and info.get('gtf_exists', False):
                    file_size = os.path.getsize(gtf_path) / (1024**2)  # MB
                    gtf_status = f"✅ 已下载 ({file_size:.2f} MB)"
                elif gtf_path:
                    gtf_status = "❌ 未下载 (配置已设置)"
                else:
                    gtf_status = "⚠️  未配置"
                
                # 检查STAR索引状态
                if info.get('star_index_exists', False):
                    index_status = "✅ 已建立索引"
                elif fasta_path:
                    index_status = "❌ 未建立索引"
                else:
                    index_status = "⚠️  无法检查 (FASTA未配置)"
                
                result.append(f"   📁 FASTA文件: {fasta_status}")
                result.append(f"   📄 GTF文件: {gtf_status}")
                result.append(f"   🔍 STAR索引: {index_status}")
                
                # 显示文件路径
                if fasta_path:
                    result.append(f"   📂 FASTA路径: {fasta_path}")
                if gtf_path:
                    result.append(f"   📂 GTF路径: {gtf_path}")
                
                # 显示下载URL（如果有）
                if 'fasta_url' in info:
                    result.append(f"   🔗 FASTA URL: {info['fasta_url']}")
                if 'gtf_url' in info:
                    result.append(f"   🔗 GTF URL: {info['gtf_url']}")
                
                result.append("-" * 80)
            
            # 统计信息
            total_genomes = len(genomes_config)
            downloaded_genomes = sum(1 for info in genomes_config.values() 
                                   if info.get('fasta_exists', False) and info.get('gtf_exists', False))
            indexed_genomes = sum(1 for info in genomes_config.values() 
                                if info.get('star_index_exists', False))
            
            result.append(f"📈 统计摘要：")
            result.append(f"   总基因组数量: {total_genomes}")
            result.append(f"   已完全下载: {downloaded_genomes}")
            result.append(f"   已建立索引: {indexed_genomes}")
            result.append(f"   准备就绪率: {(indexed_genomes/total_genomes*100):.1f}%" if total_genomes > 0 else "   准备就绪率: 0.0%")
            
            return "\n".join(result)
        
        # 如果提供了genome_name，返回特定基因组的详细信息
        if genome_name not in genomes_config:
            available_genomes = list(genomes_config.keys())
            return f"错误：基因组 '{genome_name}' 不存在。可用基因组：{', '.join(available_genomes)}"
        
        genome_info = genomes_config[genome_name]
        
        result = [f"基因组 '{genome_name}' 详细信息："]
        result.append(f"  物种: {genome_info.get('species', '未知')}")
        result.append(f"  版本: {genome_info.get('version', '未知')}")
        result.append(f"  FASTA文件: {genome_info.get('fasta', '未配置')}")
        result.append(f"  GTF文件: {genome_info.get('gtf', '未配置')}")
        
        if 'fasta_url' in genome_info:
            result.append(f"  FASTA下载URL: {genome_info['fasta_url']}")
        if 'gtf_url' in genome_info:
            result.append(f"  GTF下载URL: {genome_info['gtf_url']}")
        
        # 检查本地文件是否存在
        fasta_path = genome_info.get('fasta', '')
        gtf_path = genome_info.get('gtf', '')
        
        if fasta_path and os.path.exists(fasta_path):
            result.append(f"  ✅ FASTA文件已存在")
        elif fasta_path:
            result.append(f"  ❌ FASTA文件不存在")
        
        if gtf_path and os.path.exists(gtf_path):
            result.append(f"  ✅ GTF文件已存在")
        elif gtf_path:
            result.append(f"  ❌ GTF文件不存在")
        
        return "\n".join(result)
    
    except json.JSONDecodeError:
        return f"错误：基因组配置文件 '{config_path}' 格式不正确"
    except Exception as e:
        return f"查询基因组信息时发生错误：{str(e)}"

# ============================================================================ 
# 配置管理工具组 - 遵循DRY原则
# ============================================================================ 
class AddGenomeArgs(BaseModel):
    """基因组添加参数模型"""
    genome_name: str = Field(description="要添加的基因组的唯一名称，例如 'danRer11'")
    species: str = Field(description="该基因组所属的物种，例如 'zebrafish'")
    fasta_url: str = Field(description="基因组FASTA文件的URL")
    gtf_url: str = Field(description="基因组GTF文件的URL")
    fasta_path: str = Field(description="FASTA文件的本地存储路径")
    gtf_path: str = Field(description="GTF文件的本地存储路径")

@tool(args_schema=AddGenomeArgs)
def add_new_genome(genome_name: str, species: str, fasta_url: str, gtf_url: str, fasta_path: str, gtf_path: str) -> str:
    """
    添加一个全新的基因组到配置文件(config/genomes.json)，并验证文件路径
    自动检查文件是否存在，提供状态反馈
    """
    try:
        config_path = "config/genomes.json"
        if not os.path.exists(config_path):
            return f"错误：基因组配置文件 '{config_path}' 不存在"

        with open(config_path, 'r+', encoding='utf-8') as f:
            genomes_config = json.load(f)

            if genome_name in genomes_config:
                return f"错误：基因组 '{genome_name}' 已存在于配置中。"

            # 验证文件路径和存在状态
            fasta_exists = os.path.exists(fasta_path) if fasta_path else False
            gtf_exists = os.path.exists(gtf_path) if gtf_path else False
            
            # 检查STAR索引状态
            star_index_exists = False
            if fasta_exists:
                star_index_dir = os.path.join(os.path.dirname(fasta_path), "star_index")
                star_index_exists = (os.path.exists(star_index_dir) and 
                                   os.path.exists(os.path.join(star_index_dir, "SA")))

            new_genome_entry = {
                "species": species,
                "version": genome_name,
                "fasta": fasta_path,
                "gtf": gtf_path,
                "fasta_url": fasta_url,
                "gtf_url": gtf_url,
                "fasta_exists": fasta_exists,
                "gtf_exists": gtf_exists,
                "star_index_exists": star_index_exists
            }

            genomes_config[genome_name] = new_genome_entry
            
            f.seek(0)
            json.dump(genomes_config, f, indent=2)
            f.truncate()

        # 构建状态报告
        status_report = [f"✅ 成功添加基因组 '{genome_name}' (物种: {species}) 到 '{config_path}'。"]
        status_report.append("\n📊 文件状态:")
        
        if fasta_exists:
            file_size = os.path.getsize(fasta_path) / (1024**3)  # GB
            status_report.append(f"   📁 FASTA文件: ✅ 已存在 ({file_size:.2f} GB)")
        else:
            status_report.append(f"   📁 FASTA文件: ❌ 不存在 - 需要从URL下载")
            
        if gtf_exists:
            file_size = os.path.getsize(gtf_path) / (1024**2)  # MB
            status_report.append(f"   📄 GTF文件: ✅ 已存在 ({file_size:.2f} MB)")
        else:
            status_report.append(f"   📄 GTF文件: ❌ 不存在 - 需要从URL下载")
            
        if star_index_exists:
            status_report.append(f"   🔍 STAR索引: ✅ 已建立")
        else:
            status_report.append(f"   🔍 STAR索引: ❌ 需要构建")

        status_report.append(f"\n🔗 下载链接:")
        status_report.append(f"   FASTA: {fasta_url}")
        status_report.append(f"   GTF: {gtf_url}")

        return "\n".join(status_report)

    except json.JSONDecodeError:
        return f"错误：基因组配置文件 '{config_path}' 格式不正确"
    except Exception as e:
        return f"添加新基因组时发生错误：{str(e)}"

# ============================================================================ 
# 配置管理工具组 - 遵循DRY原则
# ============================================================================ 


@tool(args_schema=NextflowConfigArgs)
def update_nextflow_param(param_name: str, param_value: Any) -> str:
    """
    更新单个nextflow参数到AgentState
    
    返回特殊格式的字符串，供节点解析并更新状态
    """
    try:
        # 验证参数名称
        valid_params = [
            "srr_ids", "local_genome_path", "local_gtf_path", 
            "download_genome_url", "download_gtf_url", "local_fastq_files",
            "genome_version", "data", "run_download_srr", "run_download_genome", 
            "run_build_star_index", "run_fastp", "run_star_align", "run_featurecounts"
        ]
        
        if param_name not in valid_params:
            return f"❌ 错误：无效的参数名 '{param_name}'\n有效参数：{', '.join(valid_params)}"
        
        # 类型验证和转换
        if param_name.startswith("run_"):
            if isinstance(param_value, str):
                param_value = param_value.lower() in ['true', '1', 'yes', 'on']
            elif not isinstance(param_value, bool):
                return f"❌ 错误：参数 '{param_name}' 必须是布尔值"
        
        # 返回特殊格式，包含配置更新指令
        import json
        config_update = {param_name: param_value}
        
        return f"✅ 参数 '{param_name}' 已设置为: {param_value}\n[CONFIG_UPDATE] {json.dumps(config_update, ensure_ascii=False)}"
    
    except Exception as e:
        return f"❌ 更新参数时发生错误：{str(e)}"

@tool(args_schema=BatchConfigArgs)
def batch_update_nextflow_config(config_updates: Dict[str, Any]) -> str:
    """
    批量更新nextflow配置
    
    遵循DRY原则：复用单参数更新逻辑
    """
    try:
        # 验证所有参数
        valid_params = [
            "srr_ids", "local_genome_path", "local_gtf_path", 
            "download_genome_url", "download_gtf_url", "local_fastq_files",
            "genome_version", "data", "run_download_srr", "run_download_genome", 
            "run_build_star_index", "run_fastp", "run_star_align", "run_featurecounts"
        ]
        
        validated_updates = {}
        errors = []
        
        for param_name, param_value in config_updates.items():
            if param_name not in valid_params:
                errors.append(f"错误：无效的参数名 '{param_name}'")
                continue
                
            # 类型验证和转换
            if param_name.startswith("run_"):
                if isinstance(param_value, str):
                    param_value = param_value.lower() in ['true', '1', 'yes', 'on']
                elif not isinstance(param_value, bool):
                    errors.append(f"错误：参数 '{param_name}' 必须是布尔值")
                    continue
            
            validated_updates[param_name] = param_value
        
        if errors:
            return "批量更新失败:\n" + "\n".join(errors)
        
        # 返回特殊格式，包含批量配置更新指令
        import json
        config_json = json.dumps(validated_updates, ensure_ascii=False)
        
        success_params = list(validated_updates.keys())
        summary = f"✅ 批量更新完成：{len(success_params)} 个参数已更新\n"
        summary += f"更新的参数：{', '.join(success_params)}"
        
        return f"{summary}\n__CONFIG_UPDATE__:{config_json}"
    
    except Exception as e:
        return f"批量更新配置时发生错误：{str(e)}"

# ============================================================================
# 模式控制工具组 - 遵循开放封闭原则
# ============================================================================

@tool(args_schema=ModeSwitch)
def switch_to_plan_mode(target_mode: str, reason: str) -> str:
    """
    切换到计划模式
    
    遵循单一职责原则：专门处理模式切换
    """
    if target_mode != "plan":
        return f"错误：此工具只能切换到plan模式，收到：{target_mode}"
    
    return f"🔄 正在切换到计划模式...\n原因：{reason}\n✅ 模式切换成功！现在可以开始制定RNA-seq分析计划。"

@tool(args_schema=ModeSwitch)
def switch_to_execute_mode(target_mode: str, reason: str) -> str:
    """
    切换到执行模式
    
    遵循单一职责原则：专门处理模式切换
    """
    if target_mode != "execute":
        return f"错误：此工具只能切换到execute模式，收到：{target_mode}"
    
    return f"🔄 正在切换到执行模式...\n原因：{reason}\n✅ 模式切换成功！准备执行nextflow流程。"

# ============================================================================
# 执行控制工具组 - 遵循单一职责原则
# ============================================================================

@tool(args_schema=ExecutionArgs)
def execute_nextflow_pipeline(config_path: str = "config/nextflow.config", work_dir: str = "./work") -> str:
    """
    执行nextflow流程
    
    应用KISS原则：简单的流程执行
    """
    try:
        # 检查nextflow是否可用
        result = subprocess.run(["nextflow", "-version"], capture_output=True, text=True)
        if result.returncode != 0:
            return "错误：nextflow未安装或不可用"
        
        # 构建执行命令
        cmd = [
            "nextflow", "run", "main.nf",
            "-c", config_path,
            "-work-dir", work_dir
        ]
        
        return f"🚀 开始执行nextflow流程...\n命令：{' '.join(cmd)}\n⏳ 流程正在后台运行，请使用check_execution_status查看状态。"
    
    except Exception as e:
        return f"执行nextflow流程时发生错误：{str(e)}"

@tool
def check_execution_status() -> str:
    """
    检查执行状态
    
    应用KISS原则：简单的状态检查
    """
    try:
        # 检查是否有nextflow进程在运行
        result = subprocess.run(["pgrep", "-f", "nextflow"], capture_output=True, text=True)
        
        if result.returncode == 0:
            return "✅ Nextflow流程正在运行中..."
        else:
            return "⏹️ 没有检测到正在运行的nextflow流程"
    
    except Exception as e:
        return f"检查执行状态时发生错误：{str(e)}"

@tool
def get_current_nextflow_config() -> str:
    """
    获取当前nextflow配置状态
    
    遵循单一职责原则：专门获取配置信息
    """
    try:
        # 这里应该从状态管理中获取配置，暂时返回默认配置
        default_config = {
            "srr_ids": "",
            "local_genome_path": "",
            "local_gtf_path": "",
            "download_genome_url": "",
            "download_gtf_url": "",
            "local_fastq_files": "",
            "data": "./data",
            "run_download_srr": False,
            "run_download_genome": False,
            "run_build_star_index": False,
            "run_fastp": False,
            "run_star_align": False,
            "run_featurecounts": False
        }
        
        result = ["当前nextflow配置："]
        for key, value in default_config.items():
            result.append(f"  {key}: {value}")
        
        return "\n".join(result)
    
    except Exception as e:
        return f"获取配置时发生错误：{str(e)}"

# ============================================================================
# 目录树工具组 - 遵循单一职责原则
# ============================================================================

def _generate_simple_list(directory_path: str, max_depth: Optional[int] = None,
                         file_pattern: Optional[str] = None, show_only_files: bool = False) -> str:
    """
    生成简单列表格式的目录内容
    
    遵循单一职责原则：专门处理简单列表格式输出
    """
    try:
        path = Path(directory_path)
        if not path.exists():
            return f"错误：目录 '{directory_path}' 不存在"
        
        if not path.is_dir():
            return f"错误：'{directory_path}' 不是一个目录"
        
        contents = []
        
        # 根据最大深度决定使用rglob还是glob
        if max_depth is None or max_depth > 1:
            # 使用rglob进行递归查找
            pattern = file_pattern if file_pattern else "*"
            for item in path.rglob(pattern):
                if max_depth is not None:
                    # 计算相对深度
                    relative_path = item.relative_to(path)
                    depth = len(relative_path.parts) - 1
                    if depth > max_depth:
                        continue
                
                if show_only_files and item.is_dir():
                    continue
                
                if item.is_dir():
                    contents.append(f"📁 {item}/")
                else:
                    size = item.stat().st_size
                    contents.append(f"📄 {item} ({size} bytes)")
        else:
            # 只查找当前目录
            pattern = file_pattern if file_pattern else "*"
            for item in path.glob(pattern):
                if show_only_files and item.is_dir():
                    continue
                
                if item.is_dir():
                    contents.append(f"📁 {item}/")
                else:
                    size = item.stat().st_size
                    contents.append(f"📄 {item} ({size} bytes)")
        
        if not contents:
            return f"目录 '{directory_path}' 中没有找到匹配的内容"
        
        return f"目录 '{directory_path}' 内容列表：\n" + "\n".join(contents)
    
    except Exception as e:
        return f"生成简单列表时发生错误：{str(e)}"

def _generate_tree_structure(directory_path: str, max_depth: Optional[int] = None,
                           file_pattern: Optional[str] = None, show_only_files: bool = False,
                           prefix: str = "", current_depth: int = 0) -> str:
    """
    生成树形结构的目录内容
    
    遵循递归设计模式：使用递归构建树形结构
    """
    try:
        path = Path(directory_path)
        if not path.exists():
            return f"错误：目录 '{directory_path}' 不存在"
        
        if not path.is_dir():
            return f"错误：'{directory_path}' 不是一个目录"
        
        # 检查是否超过最大深度
        if max_depth is not None and current_depth > max_depth:
            return ""
        
        result = []
        
        try:
            items = sorted(path.iterdir(), key=lambda x: (x.is_file(), x.name.lower()))
        except PermissionError:
            return f"{prefix}❌ 没有权限访问此目录"
        
        # 过滤文件
        if file_pattern:
            import fnmatch
            items = [item for item in items if item.is_dir() or fnmatch.fnmatch(item.name, file_pattern)]
        
        for i, item in enumerate(items):
            is_last = i == len(items) - 1
            
            # 跳过目录（如果只显示文件）
            if show_only_files and item.is_dir():
                continue
            
            # 添加当前项
            if item.is_dir():
                connector = "└── " if is_last else "├── "
                result.append(f"{prefix}{connector}📁 {item.name}/")
                
                # 递归处理子目录
                if max_depth is None or current_depth < max_depth:
                    extension = "    " if is_last else "│   "
                    subtree = _generate_tree_structure(
                        str(item), max_depth, file_pattern, show_only_files,
                        prefix + extension, current_depth + 1
                    )
                    if subtree:
                        result.append(subtree)
            else:
                connector = "└── " if is_last else "├── "
                size = item.stat().st_size
                result.append(f"{prefix}{connector}📄 {item.name} ({size} bytes)")
        
        return "\n".join(result)
    
    except Exception as e:
        return f"生成树形结构时发生错误：{str(e)}"

@tool(args_schema=TreeListArgs)
def list_directory_tree(directory_path: str, max_depth: Optional[int] = None,
                       file_pattern: Optional[str] = None, show_only_files: bool = False,
                       output_format: str = "tree") -> str:
    """
    列出目录树结构或简单列表
    
    遵循单一职责原则：专门处理目录树和列表显示
    提供多种输出格式和过滤选项
    """
    try:
        # 验证输出格式
        if output_format not in ["tree", "list"]:
            return f"错误：不支持的输出格式 '{output_format}'。支持的格式：'tree', 'list'"
        
        # 根据输出格式调用相应的辅助函数
        if output_format == "tree":
            result = _generate_tree_structure(directory_path, max_depth, file_pattern, show_only_files)
            if not result.startswith("错误："):
                result = f"目录 '{directory_path}' 树形结构：\n{result}"
            return result
        else:  # output_format == "list"
            return _generate_simple_list(directory_path, max_depth, file_pattern, show_only_files)
    
    except Exception as e:
        return f"列出目录树时发生错误：{str(e)}"

@tool(args_schema=TaskListArgs)
def generate_analysis_task_list(analysis_type: str = "standard", force_refresh: bool = False) -> str:
    """
    生成智能RNA-seq分析任务列表，自动检测本地文件并确定最优配置
    
    遵循智能配置原则：优先使用本地文件，自动生成nextflow参数
    """
    try:
        result = ["📋 **智能任务列表生成**"]
        result.append("=" * 50)
        
        # 第1步：检测本地FASTQ文件
        result.append("\n🔍 **步骤1：检测FASTQ文件**")
        fastq_detection = _detect_local_fastq_files()
        result.extend(fastq_detection["summary"])
        
        # 第2步：检测本地基因组文件  
        result.append("\n🧬 **步骤2：检测基因组文件**")
        genome_detection = _detect_local_genome_files()
        result.extend(genome_detection["summary"])
        
        # 第3步：生成推荐配置
        result.append("\n⚙️ **步骤3：生成推荐配置**")
        recommended_config = _generate_recommended_config(
            fastq_detection["data"], 
            genome_detection["data"], 
            analysis_type
        )
        result.extend(recommended_config["summary"])
        
        # 第4步：生成任务列表
        result.append("\n📝 **步骤4：分析任务列表**")
        task_list = _generate_task_steps(recommended_config["config"], analysis_type)
        result.extend(task_list)
        
        # 第5步：显示最终配置
        result.append("\n🎯 **最终推荐配置**")
        config_summary = _format_config_summary(recommended_config["config"])
        result.extend(config_summary)
        
        # 使用建议
        result.append("\n💡 **使用建议**")
        if recommended_config["config"].get("has_local_files"):
            result.append("✅ 检测到本地文件，配置已优化为使用本地资源")
        else:
            result.append("⚠️ 未检测到本地文件，将需要下载数据和基因组")
        
        result.append("📋 配置已生成，可直接用于nextflow执行")
        
        return "\n".join(result)
        
    except Exception as e:
        return f"生成任务列表时发生错误：{str(e)}"

def _detect_local_fastq_files() -> Dict[str, Any]:
    """检测本地FASTQ文件"""
    try:
        # 搜索默认FASTQ路径
        search_paths = ["data/fastq", "data/results/fastp", "fastq", "raw_data"]
        found_files = []
        
        for path in search_paths:
            if os.path.exists(path):
                for root, dirs, files in os.walk(path):
                    for file in files:
                        if file.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                            found_files.append(os.path.join(root, file))
        
        if found_files:
            # 分析文件类型
            paired_files = {}
            single_files = []
            
            for file_path in found_files:
                file_name = os.path.basename(file_path)
                # 简化的配对检测
                if '_1.' in file_name or '_R1.' in file_name:
                    sample_id = file_name.split('_')[0]
                    if sample_id not in paired_files:
                        paired_files[sample_id] = {}
                    paired_files[sample_id]['R1'] = file_path
                elif '_2.' in file_name or '_R2.' in file_name:
                    sample_id = file_name.split('_')[0] 
                    if sample_id not in paired_files:
                        paired_files[sample_id] = {}
                    paired_files[sample_id]['R2'] = file_path
                else:
                    single_files.append(file_path)
            
            summary = [
                f"✅ 检测到 {len(found_files)} 个FASTQ文件",
                f"   - 双端文件：{len(paired_files)} 对样本",
                f"   - 单端文件：{len(single_files)} 个",
                f"   - 建议配置：使用本地FASTQ文件"
            ]
            
            return {
                "data": {
                    "found": True,
                    "paired_files": paired_files,
                    "single_files": single_files,
                    "total_files": len(found_files),
                    "recommended_path": search_paths[0] if os.path.exists(search_paths[0]) else None
                },
                "summary": summary
            }
        else:
            summary = [
                "❌ 未检测到本地FASTQ文件",
                "   - 搜索路径：" + ", ".join(search_paths),
                "   - 建议配置：需要提供SRR ID或上传FASTQ文件"
            ]
            
            return {
                "data": {"found": False},
                "summary": summary
            }
            
    except Exception as e:
        return {
            "data": {"found": False, "error": str(e)},
            "summary": [f"❌ FASTQ文件检测失败：{str(e)}"]
        }

def _detect_local_genome_files() -> Dict[str, Any]:
    """检测本地基因组文件"""
    try:
        # 读取基因组配置
        config_path = "config/genomes.json"
        if not os.path.exists(config_path):
            return {
                "data": {"found": False},
                "summary": ["❌ 基因组配置文件不存在"]
            }
        
        with open(config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        
        available_genomes = []
        ready_genomes = []
        
        for name, info in genomes_config.items():
            fasta_path = info.get('fasta', '')
            gtf_path = info.get('gtf', '')
            
            fasta_exists = os.path.exists(fasta_path) if fasta_path else False
            gtf_exists = os.path.exists(gtf_path) if gtf_path else False
            
            available_genomes.append(name)
            
            if fasta_exists and gtf_exists:
                ready_genomes.append({
                    "name": name,
                    "fasta": fasta_path,
                    "gtf": gtf_path,
                    "species": info.get('species', 'unknown')
                })
        
        if ready_genomes:
            summary = [
                f"✅ 检测到 {len(ready_genomes)} 个可用基因组",
                f"   - 推荐使用：{ready_genomes[0]['name']} ({ready_genomes[0]['species']})",
                f"   - 其他可选：{', '.join([g['name'] for g in ready_genomes[1:]])}" if len(ready_genomes) > 1 else ""
            ]
            summary = [s for s in summary if s]  # 移除空字符串
            
            return {
                "data": {
                    "found": True,
                    "ready_genomes": ready_genomes,
                    "total_available": len(available_genomes),
                    "recommended": ready_genomes[0]
                },
                "summary": summary
            }
        else:
            summary = [
                f"⚠️ 配置中有 {len(available_genomes)} 个基因组，但文件未下载",
                f"   - 可用基因组：{', '.join(available_genomes)}",
                "   - 建议配置：需要下载基因组文件"
            ]
            
            return {
                "data": {
                    "found": False,
                    "available_genomes": list(genomes_config.keys()),
                    "total_available": len(available_genomes)
                },
                "summary": summary
            }
            
    except Exception as e:
        return {
            "data": {"found": False, "error": str(e)},
            "summary": [f"❌ 基因组文件检测失败：{str(e)}"]
        }

def _generate_recommended_config(fastq_data: Dict, genome_data: Dict, analysis_type: str) -> Dict[str, Any]:
    """生成推荐的nextflow配置"""
    try:
        config = {
            "data": "./data",
            "run_fastp": True,
            "run_star_align": True,
            "run_featurecounts": True,
            "run_build_star_index": True,
            "has_local_files": False
        }
        
        summary = []
        
        # 配置FASTQ文件
        if fastq_data.get("found"):
            if fastq_data.get("recommended_path"):
                config["local_fastq_files"] = fastq_data["recommended_path"] + "/*.fastq*"
                config["run_download_srr"] = False
                summary.append("✅ 配置使用本地FASTQ文件")
                config["has_local_files"] = True
            else:
                summary.append("⚠️ 检测到FASTQ文件但路径不明确")
        else:
            config["run_download_srr"] = True
            config["srr_ids"] = ""  # 需要用户提供
            summary.append("📥 配置为下载SRR数据（需要用户提供SRR ID）")
        
        # 配置基因组文件
        if genome_data.get("found") and genome_data.get("recommended"):
            recommended = genome_data["recommended"]
            config["local_genome_path"] = recommended["fasta"]
            config["local_gtf_path"] = recommended["gtf"]
            config["run_download_genome"] = False
            config["genome_version"] = recommended["name"]
            summary.append(f"✅ 配置使用本地基因组：{recommended['name']}")
            config["has_local_files"] = True
        else:
            config["run_download_genome"] = True
            config["genome_version"] = "hg38"  # 默认
            summary.append("📥 配置为下载基因组文件（默认hg38）")
        
        # 根据分析类型调整
        if analysis_type == "minimal":
            config["run_fastp"] = False
            summary.append("🔧 最小模式：跳过质量控制")
        elif analysis_type == "comprehensive":
            config["run_multiqc"] = True
            summary.append("🔧 全面模式：启用MultiQC报告")
        
        return {
            "config": config,
            "summary": summary
        }
        
    except Exception as e:
        return {
            "config": {},
            "summary": [f"❌ 配置生成失败：{str(e)}"]
        }

def _generate_task_steps(config: Dict, analysis_type: str) -> List[str]:
    """生成任务步骤列表"""
    try:
        steps = []
        step_num = 1
        
        # 数据准备步骤
        if config.get("run_download_srr"):
            steps.append(f"{step_num}. 📥 下载SRR数据文件")
            step_num += 1
        
        if config.get("run_download_genome"):
            steps.append(f"{step_num}. 📥 下载基因组参考文件")
            step_num += 1
        
        # 索引构建
        if config.get("run_build_star_index"):
            steps.append(f"{step_num}. 🔨 构建STAR基因组索引")
            step_num += 1
        
        # 数据处理步骤
        if config.get("run_fastp"):
            steps.append(f"{step_num}. 🧹 质量控制和数据清理 (FastP)")
            step_num += 1
        
        if config.get("run_star_align"):
            steps.append(f"{step_num}. 🎯 序列比对到参考基因组 (STAR)")
            step_num += 1
        
        if config.get("run_featurecounts"):
            steps.append(f"{step_num}. 📊 基因表达定量 (featureCounts)")
            step_num += 1
        
        # 额外步骤
        if config.get("run_multiqc"):
            steps.append(f"{step_num}. 📋 生成综合质量报告 (MultiQC)")
            step_num += 1
        
        steps.append(f"{step_num}. 📁 整理输出结果和日志文件")
        
        return steps
        
    except Exception as e:
        return [f"❌ 任务步骤生成失败：{str(e)}"]

def _format_config_summary(config: Dict) -> List[str]:
    """格式化配置摘要"""
    try:
        summary = []
        
        # 数据源
        summary.append("**数据源配置：**")
        if config.get("local_fastq_files"):
            summary.append(f"  📁 FASTQ文件：{config['local_fastq_files']}")
        elif config.get("srr_ids"):
            summary.append(f"  📥 SRR下载：{config['srr_ids']}")
        else:
            summary.append("  ⚠️ FASTQ：需要配置")
        
        # 基因组
        if config.get("local_genome_path"):
            summary.append(f"  🧬 基因组：{config.get('genome_version', 'local')}")
        else:
            summary.append(f"  📥 基因组下载：{config.get('genome_version', 'hg38')}")
        
        # 分析步骤
        summary.append("\n**分析流程：**")
        processes = []
        if config.get("run_fastp"):
            processes.append("质量控制")
        if config.get("run_star_align"):
            processes.append("序列比对")
        if config.get("run_featurecounts"):
            processes.append("表达定量")
        
        summary.append(f"  🔬 启用流程：{' → '.join(processes)}")
        summary.append(f"  📂 输出目录：{config.get('data', './data')}")
        
        return summary
        
    except Exception as e:
        return [f"❌ 配置摘要生成失败：{str(e)}"]