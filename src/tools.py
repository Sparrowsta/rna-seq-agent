"""
RNA-seq智能分析助手工具模块
提供FASTQ文件查询、基因组管理、用户意图收集等核心功能
"""

import os
import json
import glob
import re
import time
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
        
        # 构建增强的结果格式
        total_samples = sum(len(samples) for samples in samples_by_dir.values())
        total_size_mb = sum(
            sum(file_info['size_mb'] for file_info in [files['R1'], files['R2'], files['single']] if file_info)
            for samples in samples_by_dir.values()
            for files in samples.values()
        )
        
        # 检测测序类型
        paired_count = sum(
            1 for samples in samples_by_dir.values()
            for files in samples.values()
            if files['R1'] and files['R2']
        )
        single_count = sum(
            1 for samples in samples_by_dir.values()
            for files in samples.values()
            if files['single']
        )
        
        # 生成智能概览
        result = f"📊 **FASTQ数据概览**\n\n"
        result += f"📈 **统计信息:**\n"
        result += f"   - 样本数量: {total_samples} 个\n"
        result += f"   - 数据大小: {total_size_mb:.1f} MB\n"
        result += f"   - 分布目录: {len(samples_by_dir)} 个\n"
        
        if paired_count > 0 and single_count == 0:
            result += f"   - 测序类型: 双端测序 (Paired-end)\n"
            sequencing_type = "双端测序"
        elif single_count > 0 and paired_count == 0:
            result += f"   - 测序类型: 单端测序 (Single-end)\n"
            sequencing_type = "单端测序"
        else:
            result += f"   - 测序类型: 混合类型 ({paired_count}双端 + {single_count}单端)\n"
            sequencing_type = "混合类型"
        
        # 估算读数和分析建议
        estimated_reads = total_size_mb * 4  # 粗略估算
        result += f"   - 预估读数: ~{estimated_reads:.0f}M reads\n"
        
        # 智能分析建议
        result += f"\n💡 **分析建议:**\n"
        if total_samples >= 3:
            result += f"   - 样本数量充足，适合差异表达分析\n"
        elif total_samples == 2:
            result += f"   - 最小比较组，可进行基础差异分析\n"
        else:
            result += f"   - 单样本，适合表达谱分析或质控\n"
        
        if sequencing_type == "双端测序":
            result += f"   - 双端数据质量较高，推荐标准流程\n"
        
        if estimated_reads >= 100:
            result += f"   - 测序深度充足，支持深度分析\n"
        elif estimated_reads >= 40:
            result += f"   - 测序深度适中，满足基础分析\n"
        else:
            result += f"   - 测序深度较低，注意质量控制\n"
        
        # 详细样本信息
        result += f"\n📁 **详细样本信息:**\n"
        for directory, samples in samples_by_dir.items():
            result += f"\n📂 {directory}:\n"
            for sample_name, files in samples.items():
                if files["R1"] and files["R2"]:
                    total_sample_size = files['R1']['size_mb'] + files['R2']['size_mb']
                    result += f"   ✅ {sample_name}: 双端配对 ({total_sample_size:.1f}MB)\n"
                    result += f"      └─ R1: {files['R1']['filename']} ({files['R1']['size_mb']}MB)\n"
                    result += f"      └─ R2: {files['R2']['filename']} ({files['R2']['size_mb']}MB)\n"
                elif files["single"]:
                    result += f"   📄 {sample_name}: 单端文件 ({files['single']['size_mb']}MB)\n"
                    result += f"      └─ {files['single']['filename']}\n"
        
        # 如果找到处理过的文件，添加提示
        if processed_files:
            result += f"\n🔄 **已处理文件:** 发现 {len(processed_files)} 个已处理文件，已自动过滤\n"
        
        # 添加后续步骤建议
        result += f"\n🚀 **后续步骤:**\n"
        result += f"   • 运行 '项目概览' 查看完整项目状态\n"
        result += f"   • 运行 '分析就绪检查' 验证系统配置\n"
        result += f"   • 使用 '/plan' 开始配置分析流程\n"
        
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

def add_genome_config(user_input: str = "") -> str:
    """使用LLM智能解析用户输入并添加基因组配置
    
    Args:
        user_input: 用户输入，包含"添加基因组"和URL信息
    
    Returns:
        添加结果的文本描述
    """
    try:
        print(f"🛠️ add_genome_config被调用，用户输入: {user_input}")
        
        if not user_input.strip():
            return "请提供基因组信息，包含FASTA和GTF文件的URL"
        
        # 使用LLM解析用户输入
        from src.core import get_shared_llm
        
        llm = get_shared_llm()
        
        system_prompt = """你是基因组信息解析专家。用户会提供包含FASTA和GTF文件URL的文本，你需要提取信息并返回JSON格式：

{
  "genome_id": "基因组标识符(如hg38,mm39,ce11)",
  "species": "物种名称(如human,mouse,caenorhabditis_elegans)", 
  "version": "版本号(通常与genome_id相同)",
  "fasta_url": "FASTA文件完整URL",
  "gtf_url": "GTF文件完整URL"
}

要求：
1. 从URL中智能识别基因组版本和物种
2. 常见映射：hg->human, mm->mouse, ce->caenorhabditis_elegans, dm->drosophila, rn->rat
3. 如果无法确定，根据生物信息学常识合理推断
4. 只返回有效的JSON，不要其他解释

示例：
输入："添加基因组 https://path/ce11/ce11.fa.gz https://path/ce11/ce11.gtf.gz"
输出：{"genome_id":"ce11","species":"caenorhabditis_elegans","version":"ce11","fasta_url":"https://path/ce11/ce11.fa.gz","gtf_url":"https://path/ce11/ce11.gtf.gz"}"""

        # 使用简单的invoke方法
        response = llm.invoke(f"System: {system_prompt}\n\nHuman: 请解析：{user_input}")
        parsed_content = response.content.strip()
        
        print(f"🤖 LLM解析结果: {parsed_content}")
        
        # 提取JSON内容
        try:
            # 尝试直接解析
            genome_info = json.loads(parsed_content)
        except json.JSONDecodeError:
            # 如果失败，尝试提取JSON部分
            json_match = re.search(r'\{.*\}', parsed_content, re.DOTALL)
            if not json_match:
                return f"LLM解析失败，无法提取有效的JSON格式。原始响应：{parsed_content}"
            
            try:
                genome_info = json.loads(json_match.group())
            except json.JSONDecodeError as e:
                return f"JSON解析错误：{str(e)}\n原始响应：{parsed_content}"
        
        # 验证必需字段
        required_fields = ['genome_id', 'species', 'version', 'fasta_url', 'gtf_url']
        missing_fields = [field for field in required_fields if not genome_info.get(field)]
        if missing_fields:
            return f"LLM解析结果缺少必需字段：{missing_fields}\n解析结果：{genome_info}"
        
        print(f"✅ 解析成功: {genome_info}")
        
        # 执行实际的基因组添加逻辑
        genome_id = genome_info['genome_id']
        species = genome_info['species']
        version = genome_info['version']
        fasta_url = genome_info['fasta_url']
        gtf_url = genome_info['gtf_url']
        
        # 验证URL格式
        if not (fasta_url.startswith('http://') or fasta_url.startswith('https://')):
            return f"FASTA URL格式无效：{fasta_url}"
        if not (gtf_url.startswith('http://') or gtf_url.startswith('https://')):
            return f"GTF URL格式无效：{gtf_url}"
        
        # 生成本地路径
        fasta_path = f"data/genomes/{species}/{version}/{version}.fa"
        gtf_path = f"data/genomes/{species}/{version}/{version}.gtf"
        
        # 构建新的基因组配置
        new_genome_config = {
            "species": species,
            "version": version,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path,
            "fasta_url": fasta_url,
            "gtf_url": gtf_url
        }
        
        # 读取现有配置
        genomes_file = Path("config/genomes.json")
        if genomes_file.exists():
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
        else:
            genomes_data = {}
        
        # 检查是否已存在
        if genome_id in genomes_data:
            existing = genomes_data[genome_id]
            return f"""基因组 {genome_id} 已存在

现有配置：
- 物种：{existing.get('species', '未知')}
- 版本：{existing.get('version', '未知')}
- FASTA URL：{existing.get('fasta_url', '未设置')}
- GTF URL：{existing.get('gtf_url', '未设置')}

如需更新，请使用不同的genome_id或先删除现有配置"""
        
        # 添加新配置
        genomes_data[genome_id] = new_genome_config
        
        # 保存配置
        os.makedirs("config", exist_ok=True)
        with open(genomes_file, 'w', encoding='utf-8') as f:
            json.dump(genomes_data, f, indent=2, ensure_ascii=False)
        
        print(f"✅ 基因组配置添加成功: {genome_id}")
        
        return f"""✅ 成功添加基因组配置：{genome_id}

📋 基因组信息：
- 基因组ID：{genome_id}
- 物种：{species}
- 版本：{version}

📁 本地路径配置：
- FASTA：{fasta_path}
- GTF：{gtf_path}

🌐 下载源：
- FASTA URL：{fasta_url}
- GTF URL：{gtf_url}

💡 提示：文件将在首次运行分析时自动下载到指定路径"""
        
    except Exception as e:
        print(f"❌ add_genome_config出错: {str(e)}")
        return f"智能解析并添加基因组时出错：{str(e)}"

# ==================== 项目信息中心工具 ====================

def get_project_overview(query: str = "") -> str:
    """项目全貌概览 - 整合所有关键信息的智能仪表板"""
    try:
        result = "🎯 **项目概览仪表板**\n\n"
        
        # 1. 数据状态总览
        result += "📊 **数据状态:**\n"
        
        # 扫描FASTQ文件
        project_root = Path(".")
        fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        all_fastq_files = []
        for ext in fastq_extensions:
            all_fastq_files.extend(project_root.rglob(ext))
        
        # 过滤原始文件
        excluded_dirs = ["work", "results", "tmp"]
        processed_indicators = ["trimmed", "fastp", "cutadapt", "filtered", "processed", "qc"]
        raw_fastq_files = []
        
        for file_path in all_fastq_files:
            if not file_path.exists() or any(excluded_dir in file_path.parts for excluded_dir in excluded_dirs):
                continue
            filename_lower = file_path.name.lower()
            if not any(indicator in filename_lower for indicator in processed_indicators):
                raw_fastq_files.append(file_path)
        
        # 统计样本信息
        total_samples = 0
        total_size_mb = 0
        sequencing_type = "未检测到"
        
        if raw_fastq_files:
            # 简单样本计数和大小统计
            sample_names = set()
            for file_path in raw_fastq_files:
                filename = file_path.name
                total_size_mb += file_path.stat().st_size / 1024 / 1024
                
                if "_1." in filename or "_R1" in filename:
                    sample_name = filename.split("_1.")[0].split("_R1")[0]
                    sample_names.add(sample_name)
                    sequencing_type = "双端测序 (Paired-end)"
                elif "_2." in filename or "_R2" in filename:
                    sample_name = filename.split("_2.")[0].split("_R2")[0]
                    sample_names.add(sample_name)
                else:
                    sample_name = filename.split(".")[0]
                    sample_names.add(sample_name)
                    sequencing_type = "单端测序 (Single-end)"
            
            total_samples = len(sample_names)
        
        result += f"   - 样本数量: {total_samples} 个\n"
        result += f"   - 数据大小: {total_size_mb:.1f} MB\n"
        result += f"   - 测序类型: {sequencing_type}\n"
        
        # 2. 基因组状态
        result += "\n🧬 **基因组状态:**\n"
        genomes_file = Path("config/genomes.json")
        ready_genomes = 0
        total_genomes = 0
        
        if genomes_file.exists():
            with open(genomes_file, 'r', encoding='utf-8') as f:
                genomes_data = json.load(f)
            total_genomes = len(genomes_data)
            
            for genome_id, info in genomes_data.items():
                fasta_path = info.get('fasta_path', '')
                gtf_path = info.get('gtf_path', '')
                if fasta_path and gtf_path and Path(fasta_path).exists() and Path(gtf_path).exists():
                    ready_genomes += 1
        
        result += f"   - 可用基因组: {total_genomes} 个\n"
        result += f"   - 就绪基因组: {ready_genomes} 个\n"
        
        # 3. 历史分析
        result += "\n📈 **历史分析:**\n"
        results_dir = Path("data/results")
        analysis_count = 0
        latest_analysis = "无"
        
        if results_dir.exists():
            analysis_dirs = [d for d in results_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
            analysis_count = len(analysis_dirs)
            
            if analysis_dirs:
                # 找最新的分析
                latest_dir = max(analysis_dirs, key=lambda x: x.stat().st_mtime)
                latest_analysis = latest_dir.name
        
        result += f"   - 历史分析: {analysis_count} 次\n"
        result += f"   - 最新分析: {latest_analysis}\n"
        
        # 4. 项目健康度评估
        result += "\n💡 **项目评估:**\n"
        
        health_score = 0
        suggestions = []
        
        if total_samples > 0:
            health_score += 40
        else:
            suggestions.append("未检测到FASTQ数据文件")
        
        if ready_genomes > 0:
            health_score += 30
        else:
            suggestions.append("需要下载基因组参考文件")
        
        if total_samples >= 3:
            health_score += 20
            suggestions.append("样本数量充足，适合差异表达分析")
        elif total_samples >= 2:
            health_score += 10
            suggestions.append("样本数量较少，考虑增加重复")
        
        if sequencing_type == "双端测序 (Paired-end)":
            health_score += 10
            suggestions.append("双端测序数据，质量较高")
        
        result += f"   - 项目健康度: {health_score}/100\n"
        
        if health_score >= 80:
            result += "   - 状态: ✅ 项目就绪，可开始分析\n"
        elif health_score >= 60:
            result += "   - 状态: ⚠️ 基本就绪，建议检查配置\n"
        else:
            result += "   - 状态: ❌ 需要完善项目配置\n"
        
        # 5. 智能建议
        if suggestions:
            result += "\n🚀 **智能建议:**\n"
            for suggestion in suggestions[:3]:  # 最多显示3个建议
                result += f"   • {suggestion}\n"
        
        result += "\n💡 输入 '/plan' 开始配置分析流程"
        
        return result.strip()
        
    except Exception as e:
        return f"生成项目概览时出错: {str(e)}"

def smart_data_detection(query: str = "") -> str:
    """智能数据检测 - 自动分析FASTQ文件配对和实验设计"""
    try:
        result = "🔍 **智能数据检测报告**\n\n"
        
        # 扫描所有FASTQ文件
        project_root = Path(".")
        fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        all_fastq_files = []
        for ext in fastq_extensions:
            all_fastq_files.extend(project_root.rglob(ext))
        
        if not all_fastq_files:
            return "❌ 未检测到任何FASTQ文件"
        
        # 过滤原始文件
        excluded_dirs = ["work", "results", "tmp"]
        processed_indicators = ["trimmed", "fastp", "cutadapt", "filtered", "processed", "qc"]
        raw_fastq_files = []
        
        for file_path in all_fastq_files:
            if not file_path.exists() or any(excluded_dir in file_path.parts for excluded_dir in excluded_dirs):
                continue
            filename_lower = file_path.name.lower()
            if not any(indicator in filename_lower for indicator in processed_indicators):
                raw_fastq_files.append(file_path)
        
        if not raw_fastq_files:
            return "❌ 未检测到原始FASTQ文件（所有文件都已被处理或在工作目录中）"
        
        # 分析样本配对
        samples = {}
        naming_patterns = []
        
        for file_path in raw_fastq_files:
            filename = file_path.name
            naming_patterns.append(filename)
            
            # 检测配对模式
            if "_1." in filename or "_R1" in filename:
                sample_name = filename.split("_1.")[0].split("_R1")[0]
                read_type = "R1"
            elif "_2." in filename or "_R2" in filename:
                sample_name = filename.split("_2.")[0].split("_R2")[0]
                read_type = "R2"
            else:
                sample_name = filename.split(".")[0]
                read_type = "single"
            
            if sample_name not in samples:
                samples[sample_name] = {"R1": None, "R2": None, "single": None, "directory": str(file_path.parent)}
            
            samples[sample_name][read_type] = {
                "filename": filename,
                "size_mb": round(file_path.stat().st_size / 1024 / 1024, 2),
                "path": str(file_path)
            }
        
        # 1. 样本配对分析
        result += "📁 **样本配对分析:**\n"
        paired_samples = 0
        single_samples = 0
        incomplete_pairs = 0
        
        for sample_name, files in samples.items():
            if files["R1"] and files["R2"]:
                paired_samples += 1
                size_diff = abs(files["R1"]["size_mb"] - files["R2"]["size_mb"])
                if size_diff > files["R1"]["size_mb"] * 0.1:  # 大小差异超过10%
                    result += f"   ⚠️ {sample_name}: 配对文件大小差异较大 ({files['R1']['size_mb']}MB vs {files['R2']['size_mb']}MB)\n"
                else:
                    result += f"   ✅ {sample_name}: 完整配对 ({files['R1']['size_mb']}MB + {files['R2']['size_mb']}MB)\n"
            elif files["single"]:
                single_samples += 1
                result += f"   📄 {sample_name}: 单端文件 ({files['single']['size_mb']}MB)\n"
            else:
                incomplete_pairs += 1
                if files["R1"]:
                    result += f"   ❌ {sample_name}: 缺少R2文件\n"
                elif files["R2"]:
                    result += f"   ❌ {sample_name}: 缺少R1文件\n"
        
        # 2. 命名规范分析
        result += f"\n📝 **命名规范分析:**\n"
        result += f"   - 双端样本: {paired_samples} 个\n"
        result += f"   - 单端样本: {single_samples} 个\n"
        result += f"   - 不完整配对: {incomplete_pairs} 个\n"
        
        # 分析命名模式
        r1_patterns = [f for f in naming_patterns if "_1." in f or "_R1" in f]
        r2_patterns = [f for f in naming_patterns if "_2." in f or "_R2" in f]
        
        if r1_patterns and r2_patterns:
            result += "   - 命名格式: 标准双端命名 (R1/R2 或 1/2)\n"
        elif single_samples > 0 and paired_samples == 0:
            result += "   - 命名格式: 单端测序命名\n"
        else:
            result += "   - 命名格式: 混合或非标准命名\n"
        
        # 3. 实验设计推测
        result += f"\n🧪 **实验设计推测:**\n"
        total_samples = len(samples)
        
        if total_samples >= 6:
            result += "   - 样本规模: 大型研究 (≥6样本)\n"
            result += "   - 建议分析: 差异表达 + 功能富集 + 共表达网络\n"
        elif total_samples >= 3:
            result += "   - 样本规模: 标准研究 (3-5样本)\n"
            result += "   - 建议分析: 差异表达分析\n"
        elif total_samples == 2:
            result += "   - 样本规模: 最小比较 (2样本)\n"
            result += "   - 建议分析: 基础差异表达（统计功效有限）\n"
        else:
            result += "   - 样本规模: 单样本\n"
            result += "   - 建议分析: 表达谱分析或质控检查\n"
        
        # 4. 质量预检
        result += f"\n📊 **数据质量预检:**\n"
        
        # 文件大小分析
        all_sizes = []
        for sample_name, files in samples.items():
            if files["R1"]:
                all_sizes.append(files["R1"]["size_mb"])
            if files["R2"]:
                all_sizes.append(files["R2"]["size_mb"])
            if files["single"]:
                all_sizes.append(files["single"]["size_mb"])
        
        if all_sizes:
            avg_size = sum(all_sizes) / len(all_sizes)
            min_size = min(all_sizes)
            max_size = max(all_sizes)
            
            result += f"   - 平均文件大小: {avg_size:.1f} MB\n"
            result += f"   - 大小范围: {min_size:.1f} - {max_size:.1f} MB\n"
            
            # 大小异常检测
            if max_size / min_size > 3:  # 最大文件是最小文件的3倍以上
                result += "   ⚠️ 文件大小差异较大，建议检查数据质量\n"
            else:
                result += "   ✅ 文件大小相对均匀\n"
            
            # 根据大小估算读数
            estimated_reads = avg_size * 4  # 粗略估算：1MB ≈ 4M reads (压缩后)
            result += f"   - 预估读数: ~{estimated_reads:.1f}M reads/样本\n"
            
            if estimated_reads < 10:
                result += "   ⚠️ 读数可能较少，注意检查测序深度\n"
            elif estimated_reads > 100:
                result += "   💡 读数充足，适合深度分析\n"
            else:
                result += "   ✅ 读数适中，满足基础分析需求\n"
        
        return result.strip()
        
    except Exception as e:
        return f"智能数据检测时出错: {str(e)}"

def check_resource_readiness(query: str = "") -> str:
    """分析就绪检查 - 评估项目分析准备度"""
    try:
        result = "🔧 **分析就绪度检查**\n\n"
        
        readiness_score = 0
        max_score = 100
        issues = []
        recommendations = []
        
        # 1. 数据文件检查 (30分)
        result += "📁 **数据文件状态:**\n"
        
        project_root = Path(".")
        fastq_extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        raw_fastq_files = []
        
        for ext in fastq_extensions:
            raw_fastq_files.extend(project_root.rglob(ext))
        
        # 过滤原始文件
        excluded_dirs = ["work", "results", "tmp"]
        processed_indicators = ["trimmed", "fastp", "cutadapt", "filtered", "processed", "qc"]
        valid_fastq_files = []
        
        for file_path in raw_fastq_files:
            if not file_path.exists():
                continue
            if any(excluded_dir in file_path.parts for excluded_dir in excluded_dirs):
                continue
            filename_lower = file_path.name.lower()
            if not any(indicator in filename_lower for indicator in processed_indicators):
                valid_fastq_files.append(file_path)
        
        if valid_fastq_files:
            result += f"   ✅ 检测到 {len(valid_fastq_files)} 个FASTQ文件\n"
            readiness_score += 30
        else:
            result += "   ❌ 未检测到有效的FASTQ文件\n"
            issues.append("缺少输入数据文件")
            recommendations.append("请确保FASTQ文件存在于data目录中")
        
        # 检查文件完整性
        corrupted_files = 0
        for file_path in valid_fastq_files:
            try:
                size = file_path.stat().st_size
                if size == 0:
                    corrupted_files += 1
            except:
                corrupted_files += 1
        
        if corrupted_files > 0:
            result += f"   ⚠️ 发现 {corrupted_files} 个可能损坏的文件\n"
            issues.append(f"{corrupted_files}个文件可能损坏")
        else:
            result += "   ✅ 所有文件完整性检查通过\n"
        
        # 2. 基因组配置检查 (40分)
        result += "\n🧬 **基因组配置状态:**\n"
        
        genomes_file = Path("config/genomes.json")
        if not genomes_file.exists():
            result += "   ❌ 基因组配置文件不存在\n"
            issues.append("缺少基因组配置文件")
            recommendations.append("运行 '添加基因组' 命令配置参考基因组")
        else:
            try:
                with open(genomes_file, 'r', encoding='utf-8') as f:
                    genomes_data = json.load(f)
                
                if not genomes_data:
                    result += "   ❌ 基因组配置为空\n"
                    issues.append("基因组配置为空")
                else:
                    result += f"   ✅ 找到 {len(genomes_data)} 个已配置基因组\n"
                    readiness_score += 20
                    
                    # 检查基因组文件完整性
                    ready_genomes = []
                    for genome_id, info in genomes_data.items():
                        fasta_path = info.get('fasta_path', '')
                        gtf_path = info.get('gtf_path', '')
                        
                        fasta_exists = fasta_path and Path(fasta_path).exists()
                        gtf_exists = gtf_path and Path(gtf_path).exists()
                        
                        if fasta_exists and gtf_exists:
                            ready_genomes.append(genome_id)
                            result += f"      ✅ {genome_id}: FASTA + GTF 就绪\n"
                        else:
                            missing = []
                            if not fasta_exists:
                                missing.append("FASTA")
                            if not gtf_exists:
                                missing.append("GTF")
                            result += f"      ⚠️ {genome_id}: 缺少 {', '.join(missing)}\n"
                    
                    if ready_genomes:
                        readiness_score += 20
                        result += f"   💡 推荐使用: {', '.join(ready_genomes[:3])}\n"
                    else:
                        issues.append("没有完整的基因组文件")
                        recommendations.append("需要下载基因组FASTA和GTF文件")
                        
            except Exception as e:
                result += f"   ❌ 读取基因组配置失败: {str(e)}\n"
                issues.append("基因组配置文件损坏")
        
        # 3. 系统资源检查 (20分)
        result += "\n💻 **系统资源评估:**\n"
        
        try:
            import psutil
            
            # 内存检查
            memory = psutil.virtual_memory()
            memory_gb = memory.total / 1024**3
            available_gb = memory.available / 1024**3
            
            result += f"   - 总内存: {memory_gb:.1f} GB\n"
            result += f"   - 可用内存: {available_gb:.1f} GB\n"
            
            if available_gb >= 16:
                result += "   ✅ 内存充足，支持大型分析\n"
                readiness_score += 10
            elif available_gb >= 8:
                result += "   ⚠️ 内存适中，建议监控使用量\n"
                readiness_score += 5
                recommendations.append("监控内存使用，必要时减少并行度")
            else:
                result += "   ❌ 内存不足，可能影响分析性能\n"
                issues.append("内存不足(<8GB)")
                recommendations.append("考虑增加内存或使用较小的数据集")
            
            # 磁盘空间检查
            disk = psutil.disk_usage('.')
            disk_free_gb = disk.free / 1024**3
            
            result += f"   - 可用磁盘空间: {disk_free_gb:.1f} GB\n"
            
            if disk_free_gb >= 100:
                result += "   ✅ 磁盘空间充足\n"
                readiness_score += 10
            elif disk_free_gb >= 50:
                result += "   ⚠️ 磁盘空间适中\n"
                readiness_score += 5
                recommendations.append("监控磁盘使用，定期清理临时文件")
            else:
                result += "   ❌ 磁盘空间不足\n"
                issues.append("磁盘空间不足(<50GB)")
                recommendations.append("清理磁盘空间或使用外部存储")
                
        except ImportError:
            result += "   ⚠️ 无法检测系统资源 (psutil未安装)\n"
            recommendations.append("安装psutil库以获得更好的资源监控")
        except Exception as e:
            result += f"   ⚠️ 资源检测错误: {str(e)}\n"
        
        # 4. 配置文件检查 (10分)
        result += "\n⚙️ **配置文件状态:**\n"
        
        config_dir = Path("config")
        if config_dir.exists():
            result += "   ✅ config目录存在\n"
            readiness_score += 5
        else:
            result += "   ⚠️ config目录不存在\n"
            recommendations.append("创建config目录存放配置文件")
        
        # 检查Nextflow配置
        nextflow_config = Path("nextflow.config")
        if nextflow_config.exists():
            result += "   ✅ 发现nextflow.config文件\n"
            readiness_score += 5
        else:
            result += "   ⚠️ 未找到nextflow.config文件\n"
            recommendations.append("将在分析时自动生成Nextflow配置")
        
        # 5. 总体评估
        result += f"\n📊 **总体就绪度: {readiness_score}/{max_score} ({readiness_score}%)**\n"
        
        if readiness_score >= 80:
            result += "🟢 **状态: 完全就绪** - 可以开始分析\n"
        elif readiness_score >= 60:
            result += "🟡 **状态: 基本就绪** - 建议解决警告后开始\n"
        elif readiness_score >= 40:
            result += "🟠 **状态: 需要配置** - 请解决关键问题\n"
        else:
            result += "🔴 **状态: 未就绪** - 需要完善基础配置\n"
        
        # 问题和建议汇总
        if issues:
            result += "\n❌ **发现的问题:**\n"
            for issue in issues:
                result += f"   • {issue}\n"
        
        if recommendations:
            result += "\n💡 **改进建议:**\n"
            for rec in recommendations:
                result += f"   • {rec}\n"
        
        result += "\n🚀 准备就绪后，输入 '/plan' 开始配置分析流程"
        
        return result.strip()
        
    except Exception as e:
        return f"就绪度检查时出错: {str(e)}"

def list_analysis_history(query: str = "") -> str:
    """历史分析管理 - 浏览和管理已完成的分析"""
    try:
        result = "📈 **分析历史记录**\n\n"
        
        results_dir = Path("data/results")
        if not results_dir.exists():
            return "📭 暂无分析历史记录\n\n💡 完成首次分析后，历史记录将显示在这里"
        
        # 扫描结果目录
        analysis_dirs = []
        for item in results_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                analysis_dirs.append(item)
        
        if not analysis_dirs:
            return "📭 results目录存在但无分析记录\n\n💡 运行分析后结果将保存在data/results/目录"
        
        # 按修改时间排序（最新在前）
        analysis_dirs.sort(key=lambda x: x.stat().st_mtime, reverse=True)
        
        result += f"📂 **发现 {len(analysis_dirs)} 个分析记录:**\n\n"
        
        for i, analysis_dir in enumerate(analysis_dirs):
            if i >= 10:  # 只显示前10个最新的分析
                break
                
            dir_name = analysis_dir.name
            modification_time = analysis_dir.stat().st_mtime
            
            # 转换时间戳为可读格式
            mod_time_str = time.strftime("%Y-%m-%d %H:%M", time.localtime(modification_time))
            
            result += f"📁 **{dir_name}**\n"
            result += f"   - 分析时间: {mod_time_str}\n"
            
            # 分析目录内容
            subdirs = []
            files = []
            total_size = 0
            
            try:
                for item in analysis_dir.iterdir():
                    if item.is_dir():
                        subdirs.append(item.name)
                    else:
                        files.append(item.name)
                        try:
                            total_size += item.stat().st_size
                        except:
                            pass
                
                total_size_mb = total_size / 1024 / 1024
                result += f"   - 结果大小: {total_size_mb:.1f} MB\n"
                
                # 识别分析类型
                analysis_types = []
                if any("fastp" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("质控")
                if any("star" in subdir.lower() or "hisat" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("比对")
                if any("bam" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("比对结果")
                if any("counts" in f.lower() or "feature" in f.lower() for f in files):
                    analysis_types.append("定量")
                if any("summary" in subdir.lower() for subdir in subdirs):
                    analysis_types.append("报告")
                
                if analysis_types:
                    result += f"   - 分析步骤: {', '.join(analysis_types)}\n"
                else:
                    result += "   - 分析步骤: 未知或部分完成\n"
                
                # 检查关键结果文件
                key_files = []
                for f in files:
                    if f.endswith(('.html', '.json', '.md')):
                        key_files.append(f)
                
                if key_files:
                    result += f"   - 关键文件: {', '.join(key_files[:3])}\n"
                    if len(key_files) > 3:
                        result += f"     (还有{len(key_files)-3}个文件...)\n"
                
                # 评估分析完整性
                if "summary" in subdirs or any("report" in f.lower() for f in files):
                    result += "   - 状态: ✅ 分析完整\n"
                elif len(subdirs) >= 2:  # 至少有2个处理步骤
                    result += "   - 状态: ⚠️ 分析部分完成\n"
                else:
                    result += "   - 状态: ❌ 分析可能未完成\n"
                    
            except Exception as e:
                result += f"   - 状态: ❌ 目录访问错误: {str(e)}\n"
            
            result += "\n"
        
        # 添加统计信息
        if len(analysis_dirs) > 10:
            result += f"⏭️ 只显示了最新的10个分析，共有{len(analysis_dirs)}个历史记录\n\n"
        
        # 分析成功配置提取
        result += "🔄 **可复用配置:**\n"
        successful_configs = []
        
        for analysis_dir in analysis_dirs[:5]:  # 检查最新5个分析
            # 查找配置文件
            config_files = []
            for item in analysis_dir.rglob("*config*"):
                if item.is_file() and item.suffix in ['.json', '.yaml', '.yml']:
                    config_files.append(item)
            
            if config_files:
                successful_configs.append({
                    'name': analysis_dir.name,
                    'time': time.strftime("%m-%d %H:%M", time.localtime(analysis_dir.stat().st_mtime)),
                    'configs': len(config_files)
                })
        
        if successful_configs:
            for config in successful_configs:
                result += f"   • {config['name']} ({config['time']}) - {config['configs']}个配置文件\n"
            result += "\n💡 这些配置可以在Plan模式中复用\n"
        else:
            result += "   • 暂无可识别的配置文件\n"
        
        result += "\n🚀 输入 '/plan' 开始新的分析流程"
        
        return result.strip()
        
    except Exception as e:
        return f"获取分析历史时出错: {str(e)}"

def get_help(query: str = "") -> str:
    """获取系统帮助信息"""
    return """
🎯 RNA-seq智能分析助手 - Normal模式 (项目信息中心)

📊 **核心项目工具:**
• 项目概览 - 一键查看项目完整状态和健康度
• 历史分析记录 - 浏览已完成的分析和可复用配置

📋 **详细信息查询:**
• 查看FASTQ文件 - 详细扫描所有测序数据文件
• 查看基因组信息 - 显示可用参考基因组状态

🗄️ **基因组管理:**
• 添加基因组配置 - 智能解析并添加新的参考基因组

🚀 **开始分析:**
输入 "/plan" 进入计划模式，Plan模式将执行深度数据检测、系统资源评估并制定智能分析方案

💡 **使用建议:**
1. 首次使用建议运行 "项目概览" 了解项目状态
2. 如需详细信息，使用具体的查询工具
3. 项目了解完成后，使用 "/plan" 进入深度分析规划

🔄 **模式分工:**
- Normal模式: 快速信息查看和项目概览
- Plan模式: 深度检测、就绪评估和分析方案制定
""".strip()