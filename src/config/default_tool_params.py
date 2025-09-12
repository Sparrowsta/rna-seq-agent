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
    
    # 注意：线程数由Nextflow resource_config统一管理，不在工具参数中重复定义
    
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
    # 注意：线程数runThreadN由Nextflow resource_config统一管理，不在这里重复定义
    
    # 基础参数
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

# FeatureCounts 默认参数 - 使用 Python 风格命名，与文档要求保持一致
DEFAULT_FEATURECOUNTS_PARAMS = {
    # 注意：线程数-T由Nextflow resource_config统一管理，不在这里重复定义
    
    # 基础参数
    "-s": 0,                                   # 链特异性（0=无，1=正向，2=反向）
    "-p": False,                               # 单端测序模式（False=单端，True=双端）
    "-B": False,                               # 双端测序时不要求两端都比对
    "-C": False,                               # 不排除嵌合体reads
    
    # 特征类型和属性  
    "-t": "exon",                              # 计数特征类型
    "-g": "gene_id",                           # 用于聚合的属性
    
    # 定量策略参数
    "-M": False,                               # 是否计数多重比对reads
    "-O": False,                               # 是否允许重叠特征分配
    "--fraction": False,                       # 是否使用分数计数模式
    "-Q": 10,                                  # 最低比对质量阈值
    
    # 重叠控制参数
    "--minOverlap": 1,                         # 最小重叠碱基数
    "--fracOverlap": 0.0,                      # 最小重叠比例（0.0=不限制）
    
    # 输出控制
    "-f": False,                               # 不按特征级别输出
    "-J": False,                               # 不计数剪接位点
    
    # 高级参数（通常保持默认或根据需要调整）
    "-a": None,                                # 注释文件路径（通过Nextflow提供）
    "-F": None,                                # 注释格式（GTF/SAF，默认GTF）
    "--primary": False,                        # 仅使用主要比对
    "--ignoreDup": False,                      # 不忽略PCR重复
    "--splitOnly": False,                      # 不限制只计数分割比对
    "--nonSplitOnly": False,                   # 不限制只计数非分割比对
    "--largestOverlap": False,                 # 不使用最大重叠策略
    "--readShiftType": None,                   # reads位置偏移类型
    "--readShiftSize": None,                   # reads位置偏移大小
    "-R": None,                                # 输出read分配详情
    "--readExtension5": None,                  # 5'端延伸长度
    "--readExtension3": None,                  # 3'端延伸长度
    "--read2pos": None,                        # Read2计数位置
    "--countReadPairs": None,                  # 计数read pairs而非单个reads
    "--donotsort": None,                       # 假设输入已排序（提高性能）
    "--byReadGroup": None,                     # 按read group分别计数
    "--extraAttributes": None,                 # 输出额外的特征属性
}

# HISAT2 默认参数 - 使用 Python 风格命名
DEFAULT_HISAT2_PARAMS = {
    # 注意：线程数-p由Nextflow resource_config统一管理，不在这里重复定义
    
    # 基础RNA-seq参数
    "rna_strandness": "unstranded",            # 链特异性: unstranded/RF/FR
    "dta": True,                               # 启用下游转录本组装模式
    
    # 比对质量控制参数
    "score_min": "L,0,-0.2",                   # 最小比对得分函数
    "mp": "6,2",                               # 错配惩罚: max,min
    "sp": "2,1",                               # 软剪切惩罚: max,min
    "np": 1,                                   # N碱基惩罚
    "rdg": "5,3",                              # read gap开放和延伸惩罚
    "rfg": "5,3",                              # 参考基因组gap开放和延伸惩罚
    
    # 内含子长度控制
    "max_intronlen": 500000,                   # 最大内含子长度
    "min_intronlen": 20,                       # 最小内含子长度
    
    # 比对控制参数
    "no_mixed": False,                         # 允许混合比对（单端+成对）
    "no_discordant": False,                    # 允许不一致比对
    "k": 5,                                    # 最多报告k个比对结果
    "a": False,                                # 报告所有比对（与-k冲突）
    "n_ceil": "L,0,0.15",                      # N碱基数量上限函数
    
    # 种子和扩展参数
    "i": "S,1,1.15",                           # 种子间隔函数
    "L": 20,                                   # 种子长度
    "N": 0,                                    # 每个种子最多错配数
    "D": 15,                                   # 扩展尝试次数
    "R": 2,                                    # 重新播种次数
    
    # 输出控制参数
    "no_unal": False,                          # 不输出未比对reads
    "no_hd": False,                            # 不输出SAM头部
    "no_sq": False,                            # 不输出SQ行
    "rg_id": None,                             # Read Group ID
    "rg": None,                                # Read Group标签
    
    # 性能优化参数
    "mm": False,                               # 使用内存映射文件
    "qc_filter": False,                        # 跳过质控失败的reads
    "skip": None,                              # 跳过前N个reads
    "upto": None,                              # 最多处理N个reads
    "trim5": None,                             # 修剪5'端碱基数
    "trim3": None,                             # 修剪3'端碱基数
    
    # 高级参数（RNA-seq专用）
    "novel_splicesite_outfile": None,          # 新剪接位点输出文件
    "novel_splicesite_infile": None,           # 新剪接位点输入文件
    "no_temp_splicesite": False,               # 不使用临时剪接位点
    "no_spliced_alignment": False,             # 不进行剪接比对
    "transcriptome_mapping_only": False,       # 仅转录组比对模式
    
    # 压缩文件支持
    "phred33": None,                           # 输入质量为Phred+33
    "phred64": None,                           # 输入质量为Phred+64
    "int_quals": None,                         # 质量值为空格分隔整数
    
    # 其他输出选项
    "met_file": None,                          # 比对统计输出文件
    "met_stderr": False,                       # 比对统计输出到stderr
    "new_summary": False,                      # 使用新的统计格式
    "summary_file": None,                      # 统计摘要输出文件
    "quiet": False,                            # 安静模式
    "un": None,                                # 未比对reads输出文件前缀
    "al": None,                                # 已比对reads输出文件前缀
    "un_conc": None,                           # 未比对配对reads输出前缀
    "al_conc": None,                           # 已比对配对reads输出前缀
}
