"""
工具默认参数配置
定义 RNA-seq 流水线中各工具的默认参数
"""

# FastP 默认参数（与 fastp.nf 的 Nextflow 参数保持对齐）
DEFAULT_FASTP_PARAMS = {
    # 基础功能开关
    "adapter_trimming": True,           # 自动adapter trimming
    "quality_filtering": True,          # 质量过滤
    "length_filtering": True,           # 长度过滤
    
    # 质量参数
    "qualified_quality_phred": 20,      # 质量阈值
    "unqualified_percent_limit": 40,    # 低质量base比例限制
    "n_base_limit": 5,                  # N碱基数量限制
    
    # 长度参数
    "length_required": 15,              # 最短长度要求
    
    # 输出参数
    "html_report": True,                # 生成HTML报告
    "json_report": True,                # 生成JSON报告
    
    # 性能参数（与 Nextflow 对齐）
    # - fastp_cpus: Nextflow process.cpus；默认为4
    # - threads    : fastp 的 --thread 参数；在 fastp.nf 内部优先取 threads，否则回退 fastp_cpus，再回退 4
    "fastp_cpus": 4,
    
    # 额外 fastp 高级参数（按需传递）
    # 输入质量编码与读取控制
    "phred64": False,                   # -6/--phred64 输入为phred64质量
    "reads_to_process": None,           # --reads_to_process 仅处理前N条reads
    "fix_mgi_id": False,                # --fix_mgi_id 修复MGI测序ID
    
    # PE adapter 自动检测（Nextflow 默认启用）
    "detect_adapter_for_pe": True,      # 与 fastp.nf 的默认保持一致
    
    # 前后端定长修剪与最大长度
    "trim_front1": None,                # -f/--trim_front1
    "trim_tail1": None,                 # -t/--trim_tail1
    "max_len1": None,                   # -b/--max_len1
    "trim_front2": None,                # -F/--trim_front2
    "trim_tail2": None,                 # -T/--trim_tail2
    "max_len2": None,                   # -B/--max_len2
    
    # polyG / polyX 修剪（Nextflow 默认启用 trim_poly_g）
    "trim_poly_g": True,                # -g/--trim_poly_g（Illumina NextSeq/NovaSeq常见）
    "poly_g_min_len": None,             # --poly_g_min_len
    "disable_trim_poly_g": False,       # -G/--disable_trim_poly_g（默认关闭）
    "trim_poly_x": False,               # -x/--trim_poly_x（默认关闭）
    "poly_x_min_len": None,             # --poly_x_min_len
    
    # 滑窗切除与门限（Nextflow 默认关闭开关，窗口/门限有缺省）
    "cut_front": False,                 # -5/--cut_front（默认关闭）
    "cut_tail": False,                  # -3/--cut_tail（默认关闭）
    "cut_right": False,                 # -r/--cut_right（默认关闭）
    "cut_window_size": 4,               # -W/--cut_window_size（fastp.nf 缺省为4）
    "cut_mean_quality": 20,             # -M/--cut_mean_quality（fastp.nf 缺省为20）
    "cut_front_window_size": None,      # --cut_front_window_size
    "cut_front_mean_quality": None,     # --cut_front_mean_quality
    "cut_tail_window_size": None,       # --cut_tail_window_size
    "cut_tail_mean_quality": None,      # --cut_tail_mean_quality
    "cut_right_window_size": None,      # --cut_right_window_size
    "cut_right_mean_quality": None,     # --cut_right_mean_quality
    
    # 质量/长度过滤细化（与 fastp.nf 默认对齐）
    "average_qual": None,               # -e/--average_qual
    "disable_length_filtering": False,  # -L/--disable_length_filtering（默认关闭）
    "length_limit": None,               # --length_limit
    "low_complexity_filter": False,     # -y/--low_complexity_filter（默认关闭）
    "complexity_threshold": None,       # -Y/--complexity_threshold
    
    # PE 重叠校正与检测（默认关闭）
    "correction": False,                # -c/--correction（仅PE）
    "overlap_len_require": None,        # --overlap_len_require
    "overlap_diff_limit": None,         # --overlap_diff_limit
    "overlap_diff_percent_limit": None, # --overlap_diff_percent_limit
    
    # 过表达序列分析（Nextflow 默认开启分析，采样阈值按需）
    "overrepresentation_analysis": True, # 与 fastp.nf 默认一致
    "overrepresentation_sampling": None, # -P/--overrepresentation_sampling
}

# STAR 默认参数 - 使用 Python 风格命名
DEFAULT_STAR_PARAMS = {
    # 基础参数
    "runThreadN": 8,                          # 默认线程数
    "outSAMtype": "BAM SortedByCoordinate",   # 输出排序后的BAM
    "outSAMunmapped": "Within",               # 未比对reads保留在BAM中
    "outSAMattributes": "All",                # 输出所有SAM属性
    
    # RNA-seq 标准参数
    "outFilterMultimapNmax": 20,              # 最多允许比对到20个位置
    "alignSJoverhangMin": 8,                  # 剪接位点最小overhang
    "alignSJDBoverhangMin": 1,                # 注释剪接位点最小overhang
    "outFilterMismatchNmax": 999,             # 不限制错配数（由后续参数控制）
    "outFilterMismatchNoverReadLmax": 0.04,   # 最多4%错配率
    "alignIntronMin": 20,                     # 最小内含子长度
    "alignIntronMax": 1000000,                # 最大内含子长度
    "alignMatesGapMax": 1000000,              # 配对reads最大间隔
    
    # 定量相关参数
    "quantMode": "TranscriptomeSAM GeneCounts", # 生成转录组BAM和基因计数
    "twopassMode": "Basic",                   # 使用两遍比对模式
    
    # 内存和性能优化
    "limitBAMsortRAM": 10000000000,           # BAM排序内存限制 10GB
    "outBAMsortingThreadN": 4,                # BAM排序线程数
    "genomeLoad": "NoSharedMemory",           # 不使用共享内存
    
    # 输出文件前缀
    "outFileNamePrefix": "STAR_",             # 输出文件前缀
    
    # 额外的高级参数
    "star_cpus": None,                        # Nextflow process.cpus
    "readFilesCommand": None,                 # 读取压缩文件的命令
    "outReadsUnmapped": None,                 # 未比对reads输出选项
    "outFilterIntronMotifs": None,            # 内含子motif过滤
    "outSAMstrandField": None,                # 链信息字段
    "outFilterType": None,                    # 过滤类型
    "sjdbGTFfile": None,                      # GTF注释文件（构建索引用）
    "sjdbOverhang": None,                     # 剪接位点overhang（构建索引用）
    "chimSegmentMin": None,                   # 嵌合体最小片段长度
    "chimOutType": None,                      # 嵌合体输出类型
    "chimMainSegmentMultNmax": None,          # 嵌合体主片段最大比对数
}

# FeatureCounts 默认参数 - 使用 Python 风格命名
DEFAULT_FEATURECOUNTS_PARAMS = {
    # 基础参数
    "T": 4,                                    # 默认线程数
    "p": False,                                # 单端测序不需要配对检查
    "B": False,                                # 不要求配对reads都比对
    "C": False,                                # 不排除嵌合体reads
    
    # 特征类型和属性
    "t": "exon",                               # 计数特征类型
    "g": "gene_id",                            # 用于聚合的属性
    
    # 定量参数
    "s": 0,                                    # 非链特异性文库
    "Q": 10,                                   # 最低比对质量
    "M": True,                                 # 计数多重比对reads
    "O": True,                                 # 允许reads分配到重叠特征
    "fraction": True,                          # 多重比对reads分数计数
    
    # RNA-seq 特殊参数
    "primary": False,                          # 计数所有比对（不仅主要比对）
    "minOverlap": 10,                          # 最小重叠碱基数
    "fracOverlap": 0.0,                        # 不要求最小重叠比例
    "largestOverlap": False,                   # 不使用最大重叠分配策略
    
    # 输出格式
    "f": False,                                # 不按特征级别输出
    "J": False,                                # 不计数剪接位点
    "extraAttributes": None,                   # 不输出额外属性
    
    # 额外的高级参数
    "featurecounts_cpus": None,                # Nextflow process.cpus
    "a": None,                                 # 注释文件路径（通过Nextflow提供）
    "F": None,                                 # 注释格式（GTF/SAF）
    "R": None,                                 # 输出read分配信息
    "readExtension5": None,                    # 5'端延伸长度
    "readExtension3": None,                    # 3'端延伸长度
    "read2pos": None,                          # Read2计数位置
    "countReadPairs": None,                    # 是否计数read pairs
    "ignoreDup": None,                         # 是否忽略重复reads
    "splitOnly": None,                         # 只计数split alignment
    "nonSplitOnly": None,                      # 只计数non-split alignment
    "countSplitAlignmentsOnly": None,          # 计数split alignments
    "byReadGroup": None,                       # 按read group计数
    "donotsort": None,                         # 不排序（假定已排序）
}
