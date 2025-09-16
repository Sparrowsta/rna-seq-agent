"""
RNA-seq智能分析助手 - 工具模块

重构后的工具组织结构：
- common_tools: 通用检测和扫描工具 (scan_fastq_files, scan_system_resources, scan_genome_files)
- genome_tools: 基因组管理工具 (add_genome_config, download_genome_assets, build_star_index, build_hisat2_index)
- fastp_tools: FastP质控工具 (run_nextflow_fastp, parse_fastp_results)
- star_tools: STAR比对工具 (run_nextflow_star, parse_star_metrics)
- hisat2_tools: HISAT2比对工具 (run_nextflow_hisat2, parse_hisat2_metrics)
- featurecounts_tools: FeatureCounts定量工具 (run_nextflow_featurecounts, parse_featurecounts_metrics)
- analysis_tools: 分析报告工具 (write_analysis_markdown)
- utils_tools: 通用工具 (get_project_overview, list_analysis_history, write_params_file)
"""

# 从各个模块导入工具函数
from .common_tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files
)

from .genome_tools import (
    add_genome_config,
    download_genome_assets,
    build_star_index,
    build_hisat2_index
)

from .fastp_tools import (
    run_nextflow_fastp,
    parse_fastp_results
)

from .star_tools import (
    run_nextflow_star,
    parse_star_metrics
)

from .hisat2_tools import (
    run_nextflow_hisat2,
    parse_hisat2_metrics
)

from .featurecounts_tools import (
    run_nextflow_featurecounts,
    parse_featurecounts_metrics
)

from .analysis_tools import (
    write_analysis_markdown
)

from .utils_tools import (
    get_project_overview,
    list_analysis_history,
    write_params_file,
    enhance_tool_result_with_debug
)

# 为了兼容性，导出所有工具
__all__ = [
    # 通用检测工具
    "scan_fastq_files",
    "scan_system_resources", 
    "scan_genome_files",
    
    # 基因组管理
    "add_genome_config",
    "download_genome_assets",
    "build_star_index",
    "build_hisat2_index",
    
    # FastP质控
    "run_nextflow_fastp",
    "parse_fastp_results",
    
    # STAR比对
    "run_nextflow_star",
    "parse_star_metrics",
    
    # HISAT2比对
    "run_nextflow_hisat2",
    "parse_hisat2_metrics",
    
    # FeatureCounts定量
    "run_nextflow_featurecounts",
    "parse_featurecounts_metrics",
    
    # 分析报告
    "write_analysis_markdown",
    
    # 通用工具
    "get_project_overview",
    "list_analysis_history",
    "write_params_file",
    "enhance_tool_result_with_debug"
]

# 工具分类映射，便于节点按需导入
TOOL_CATEGORIES = {
    "detect": [
        "scan_fastq_files",
        "scan_system_resources",
        "scan_genome_files"
    ],
    "normal": [
        "scan_fastq_files",
        "scan_genome_files",
        "add_genome_config",
        "get_project_overview", 
        "list_analysis_history"
    ],
    "prepare": [
        "scan_fastq_files",
        "scan_system_resources",
        "scan_genome_files",
        "get_project_overview",
        "write_params_file"
    ],
    "fastp": [
        "run_nextflow_fastp",
        "parse_fastp_results"
    ],
    "star": [
        "download_genome_assets",
        "build_star_index",
        "run_nextflow_star", 
        "parse_star_metrics",
        "scan_genome_files"
    ],
    "hisat2": [
        "download_genome_assets",
        "build_hisat2_index",
        "run_nextflow_hisat2",
        "parse_hisat2_metrics",
        "scan_genome_files"
    ],
    "featurecounts": [
        "run_nextflow_featurecounts",
        "parse_featurecounts_metrics",
        "scan_genome_files"
    ],
    "analysis": [
        "parse_fastp_results",
        "parse_star_metrics", 
        "parse_hisat2_metrics",
        "parse_featurecounts_metrics",
        "write_analysis_markdown"
    ]
}