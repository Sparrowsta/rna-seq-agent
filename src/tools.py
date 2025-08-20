"""
RNA-seq智能分析助手工具模块
提供FASTQ文件查询、基因组管理、用户意图收集等核心功能
"""

import os
import json
import glob
from pathlib import Path
from typing import Dict, List, Any

# ==================== 信息查询工具 ====================

def query_fastq_files(query: str = "") -> str:
    """在整个工作目录递归扫描并查询FASTQ文件信息"""
    try:
        # 递归搜索当前工作目录下的所有FASTQ文件
        project_root = Path(".")
        fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        
        all_fastq_files = []
        for ext in fastq_extensions:
            all_fastq_files.extend(project_root.rglob(ext))
        
        if not all_fastq_files:
            return "在项目目录中未找到任何FASTQ文件"
        
        # 过滤掉工作目录、结果目录和处理过的文件
        excluded_dirs = ["work", "results", "tmp"]
        processed_indicators = ["trimmed", "fastp", "cutadapt", "filtered", "processed", "qc"]
        raw_fastq_files = []
        processed_files = []
        excluded_files = []
        
        for file_path in all_fastq_files:
            # 检查文件是否真实存在且可访问
            if not file_path.exists() or not file_path.is_file():
                continue
                
            # 检查是否在排除目录中
            if any(excluded_dir in file_path.parts for excluded_dir in excluded_dirs):
                excluded_files.append(file_path)
                continue
                
            filename_lower = file_path.name.lower()
            if any(indicator in filename_lower for indicator in processed_indicators):
                processed_files.append(file_path)
            else:
                raw_fastq_files.append(file_path)
        
        if not raw_fastq_files:
            return "没有找到原始FASTQ文件（所有文件都在工作目录或已被处理）"
        
        # 按目录分组原始文件
        samples_by_dir = {}
        
        for file_path in raw_fastq_files:
            directory = str(file_path.parent)
            if directory not in samples_by_dir:
                samples_by_dir[directory] = {}
            
            filename = file_path.name
            # 提取样本名称
            if "_1." in filename or "_R1" in filename:
                sample_name = filename.split("_1.")[0].split("_R1")[0]
                read_type = "R1"
            elif "_2." in filename or "_R2" in filename:
                sample_name = filename.split("_2.")[0].split("_R2")[0]
                read_type = "R2"
            else:
                sample_name = filename.split(".")[0]
                read_type = "single"
            
            if sample_name not in samples_by_dir[directory]:
                samples_by_dir[directory][sample_name] = {"R1": None, "R2": None, "single": None}
            
            samples_by_dir[directory][sample_name][read_type] = {
                "filename": filename,
                "size_mb": round(file_path.stat().st_size / 1024 / 1024, 2),
                "full_path": str(file_path)
            }
        
        # 构建结果
        total_samples = sum(len(samples) for samples in samples_by_dir.values())
        result = f"递归扫描发现 {total_samples} 个FASTQ样本（分布在 {len(samples_by_dir)} 个目录）:\n\n"
        
        for directory, samples in samples_by_dir.items():
            result += f"📂 目录: {directory}\n"
            for sample_name, files in samples.items():
                result += f"   📁 样本: {sample_name}\n"
                if files["R1"] and files["R2"]:
                    result += f"      - R1: {files['R1']['filename']} ({files['R1']['size_mb']} MB)\n"
                    result += f"      - R2: {files['R2']['filename']} ({files['R2']['size_mb']} MB)\n"
                    result += "      - 类型: 双端测序\n"
                elif files["single"]:
                    result += f"      - 文件: {files['single']['filename']} ({files['single']['size_mb']} MB)\n"
                    result += "      - 类型: 单端测序\n"
            result += "\n"
        
        # 如果找到处理过的文件，添加提示
        if processed_files:
            result += f"💡 另外发现 {len(processed_files)} 个已处理的FASTQ文件，已自动过滤\n"
        
        return result.strip()
        
    except Exception as e:
        return f"查询FASTQ文件时出错: {str(e)}"

def query_genome_info(query: str = "") -> str:
    """查询基因组配置信息"""
    try:
        genomes_file = Path("config/genomes.json")
        if not genomes_file.exists():
            return "未找到config/genomes.json配置文件"
        
        with open(genomes_file, 'r', encoding='utf-8') as f:
            genomes_data = json.load(f)
        
        if not genomes_data:
            return "基因组配置文件为空"
        
        result = f"可用基因组 ({len(genomes_data)} 个):\n\n"
        
        for genome_id, info in genomes_data.items():
            species = info.get('species', '未知物种')
            version = info.get('version', genome_id)
            
            # 检查本地文件状态
            fasta_path_str = info.get('fasta_path', '')
            gtf_path_str = info.get('gtf_path', '')
            
            # 只有路径非空且文件存在时才认为已下载
            fasta_exists = fasta_path_str and Path(fasta_path_str).exists()
            gtf_exists = gtf_path_str and Path(gtf_path_str).exists()
            
            fasta_status = "✅ 已下载" if fasta_exists else "❌ 未下载"
            gtf_status = "✅ 已下载" if gtf_exists else "❌ 未下载"
            
            result += f"🧬 {genome_id} ({species})\n"
            result += f"   - 版本: {version}\n"
            result += f"   - FASTA: {fasta_status}\n"
            result += f"   - GTF: {gtf_status}\n"
            
            if fasta_exists:
                size_mb = round(Path(fasta_path_str).stat().st_size / 1024 / 1024, 2)
                result += f"   - FASTA大小: {size_mb} MB\n"
            
            result += "\n"
        
        return result.strip()
        
    except Exception as e:
        return f"查询基因组信息时出错: {str(e)}"

def list_directory_tree(query: str = "") -> str:
    """显示data目录结构"""
    try:
        def build_tree(path: Path, prefix: str = "", max_depth: int = 3, current_depth: int = 0) -> str:
            if current_depth >= max_depth:
                return ""
            
            items = []
            try:
                for item in sorted(path.iterdir()):
                    if item.name.startswith('.'):
                        continue
                    items.append(item)
            except PermissionError:
                return f"{prefix}❌ 权限不足\n"
            
            result = ""
            for i, item in enumerate(items):
                is_last = i == len(items) - 1
                current_prefix = "└── " if is_last else "├── "
                next_prefix = prefix + ("    " if is_last else "│   ")
                
                if item.is_dir():
                    result += f"{prefix}{current_prefix}📁 {item.name}/\n"
                    if current_depth < max_depth - 1:
                        result += build_tree(item, next_prefix, max_depth, current_depth + 1)
                else:
                    # 为大文件显示大小信息
                    size_info = ""
                    try:
                        size_bytes = item.stat().st_size
                        if size_bytes > 1024 * 1024:  # >1MB
                            size_mb = round(size_bytes / 1024 / 1024, 1)
                            size_info = f" ({size_mb}M)"
                        elif size_bytes > 1024:  # >1KB
                            size_kb = round(size_bytes / 1024, 1)
                            size_info = f" ({size_kb}K)"
                    except:
                        pass
                    
                    result += f"{prefix}{current_prefix}📄 {item.name}{size_info}\n"
            
            return result
        
        data_root = Path("data")
        if not data_root.exists() or not data_root.is_dir():
            return "❌ data目录不存在"
            
        result = f"📂 data目录结构:\n\n"
        result += f"📁 data/\n"
        result += build_tree(data_root, max_depth=3)
        
        return result.strip()
        
    except Exception as e:
        return f"查看目录结构时出错: {str(e)}"

def get_help(query: str = "") -> str:
    """获取系统帮助信息"""
    # query参数保留用于工具接口一致性
    return """
🧬 RNA-seq智能分析助手 - Normal模式功能

📋 信息查询:
• 查看FASTQ文件 - 扫描所有测序数据文件
• 查看基因组信息 - 显示可用参考基因组
• 浏览目录结构 - 查看项目文件组织

🗄️ 基因组管理:
• 搜索UCSC基因组 - 从数据库查找新基因组
• 添加UCSC基因组 - 自动配置基因组信息
• 添加自定义基因组 - 手动配置基因组

🚀 开始分析:
输入 "/plan" 进入计划模式，开始配置RNA-seq分析流程
""".strip()