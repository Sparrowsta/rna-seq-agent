# agent/prompt.py

SYSTEM_PROMPT = """
你是一个生物信息学分析助手，使用React模式工作。

**React模式工作流程**:
1. **思考**: 分析当前情况，确定下一步行动
2. **行动**: 调用合适的工具执行操作  
3. **观察**: 分析工具执行结果
4. **重复**: 基于结果继续思考-行动-观察循环

**重要规则**:
- 每次决定行动前，建议先分析当前情况
- 思考是可选的，但行动是必须的
- 每次只调用一个工具
- 直接调用工具函数，不要添加额外标记
- 绝对禁止使用 `REDACTED_SPECIAL_TOKEN`、`<JSON>`、`function` 标签
- 如果下载类问题超时，自动重试3次，如果仍然失败，报告错误

**严格禁止幻觉**:
- 不能虚构工具调用结果
- 不能假设文件存在或不存在
- 所有观察结果必须基于实际工具调用
- 如果工具调用失败，必须报告真实错误信息

**智能跳过策略**:
由于生物信息学分析时间很长，进程可能中断。在每个分析步骤前都要检查文件是否存在：
- 如果文件已存在，跳过该步骤，继续下一步
- 如果文件不存在，执行该步骤
- 这样可以避免重复计算，节省时间

**检测工具**:
- `check_files_exist_tool` - 通用文件验证工具，支持检查：
  - `star_index` - 检查STAR索引是否存在
  - `gtf_file` - 检查GTF文件是否存在
  - `bam_files` - 检查BAM文件是否存在
  - `fastq_files` - 检查FASTQ文件是否存在
  - `fastp_results` - 检查FASTP结果文件是否存在
  - `counts_file` - 检查featureCounts结果文件是否存在
- `check_environment_tool` - 检查环境和工具可用性

**RNA-seq分析流程（带智能跳过）**:
1. 检查环境 (check_environment_tool)
2. 搜索基因组 (search_genome_tool)
   - 如果存在：跳过下载步骤
   - 如果不存在：下载基因组 (download_genome_tool)
3. 搜索FASTQ文件 (search_fastq_tool)
   - 如果存在：跳过下载步骤
   - 如果不存在：下载FASTQ文件 (download_fastq_tool)
4. 验证FASTQ文件 (validate_fastq_tool)
5. **检查STAR索引** (check_files_exist_tool file_type=star_index)
   - 如果存在：跳过构建步骤
   - 如果不存在：构建STAR索引 (build_star_index_tool)
6. **检查FASTP结果** (check_files_exist_tool file_type=fastp_results)
   - 如果存在：跳过质量控制
   - 如果不存在：质量控制 (run_fastp_tool)
7. **检查STAR比对结果** (check_files_exist_tool file_type=bam_files)
   - 如果存在：跳过比对
   - 如果不存在：序列比对 (run_star_align_tool)
8. **检查featureCounts结果** (check_files_exist_tool file_type=counts_file)
   - 如果存在：跳过定量
   - 如果不存在：基因定量 (run_featurecounts_tool)
9. **收集结果** (collect_results_tool)
10. **生成报告** (generate_report_tool)

**智能跳过示例**:
```
思考: 在运行STAR比对前，我需要检查BAM文件是否已经存在
行动: 调用check_files_exist_tool检查BAM文件 (file_type=bam_files)
观察: BAM文件已存在，可以跳过STAR比对步骤
思考: BAM文件已存在，直接进行下一步定量分析
行动: 调用check_files_exist_tool检查定量结果 (file_type=counts_file)
观察: 定量结果不存在，需要运行featureCounts
思考: 需要运行featureCounts进行基因定量
行动: 调用run_featurecounts_tool进行定量分析
```

**关键**: 每个观察后必须继续下一个行动循环，只有完成所有步骤才算分析完成。优先检查文件是否存在，避免重复计算。
"""

# 工具列表
TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "check_environment_tool",
            "description": "检查conda环境和工具可用性",
            "parameters": {"type": "object", "properties": {}, "required": []}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "setup_environment_tool",
            "description": "设置分析环境，安装必要工具",
            "parameters": {"type": "object", "properties": {}, "required": []}
        }
    },
    {
        "type": "function", 
        "function": {
            "name": "search_genome_tool",
            "description": "搜索和列出基因组信息 - 可查询特定基因组或列出所有可用基因组",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {"type": "string", "description": "基因组名称（可选，不提供则列出所有基因组）"}
                },
                "required": []
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "search_fastq_tool", 
            "description": "搜索和验证FASTQ文件 - 检查文件存在性并验证完整性",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_id": {"type": "string", "description": "SRR ID"}
                },
                "required": ["srr_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "download_genome_tool",
            "description": "下载基因组文件",
            "parameters": {
                "type": "object", 
                "properties": {
                    "genome_name": {"type": "string", "description": "基因组名称"}
                },
                "required": ["genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "download_fastq_tool",
            "description": "下载FASTQ文件",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_id": {"type": "string", "description": "SRR ID"}
                },
                "required": ["srr_id"]
            }
        }
    },

    {
        "type": "function",
        "function": {
            "name": "check_files_exist_tool",
            "description": "通用文件验证工具 - 检查各种类型的文件是否存在",
            "parameters": {
                "type": "object",
                "properties": {
                    "file_type": {"type": "string", "description": "文件类型: star_index, gtf_file, bam_files, fastq_files, fastp_results, counts_file"},
                    "genome_name": {"type": "string", "description": "基因组名称（用于star_index和gtf_file类型）"},
                    "srr_ids": {"type": "array", "items": {"type": "string"}, "description": "SRR ID列表（用于bam_files和fastq_files类型）"}
                },
                "required": ["file_type"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "build_star_index_tool",
            "description": "构建STAR索引",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {"type": "string", "description": "基因组名称"}
                },
                "required": ["genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_fastp_tool",
            "description": "运行质量控制",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_id": {"type": "string", "description": "SRR ID"}
                },
                "required": ["srr_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_star_align_tool",
            "description": "运行序列比对",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_id": {"type": "string", "description": "SRR ID"},
                    "genome_name": {"type": "string", "description": "基因组名称"}
                },
                "required": ["srr_id", "genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_featurecounts_tool",
            "description": "运行基因定量",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_ids": {"type": "array", "items": {"type": "string"}, "description": "SRR ID列表"},
                    "genome_name": {"type": "string", "description": "基因组名称"}
                },
                "required": ["srr_ids", "genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "collect_results_tool",
            "description": "收集分析结果并解析关键指标 - 解析fastp JSON、STAR log、featureCounts summary文件，生成详细的上游分析统计",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_ids": {"type": "array", "items": {"type": "string"}, "description": "SRR ID列表"}
                },
                "required": ["srr_ids"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "generate_report_tool",
            "description": "生成分析报告 - 包含详细的上游分析总结，自动解析fastp、STAR、featureCounts结果，生成HTML格式的综合报告",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_ids": {"type": "array", "items": {"type": "string"}, "description": "SRR ID列表"}
                },
                "required": ["srr_ids"]
            }
        }
    },


    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": "列出指定路径的文件",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {"type": "string", "description": "文件路径"},
                    "recursive": {"type": "boolean", "description": "是否递归列出子目录", "default": False}
                },
                "required": ["path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "add_genome_to_config",
            "description": "添加新基因组到配置文件",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {"type": "string", "description": "基因组名称"},
                    "species": {"type": "string", "description": "物种名称"},
                    "fasta_url": {"type": "string", "description": "FASTA文件URL"},
                    "gtf_url": {"type": "string", "description": "GTF文件URL"}
                },
                "required": ["genome_name", "species", "fasta_url", "gtf_url"]
            }
        }
    },

    {
        "type": "function",
        "function": {
            "name": "unsupported_request",
            "description": "处理不支持的请求",
            "parameters": {
                "type": "object",
                "properties": {
                    "user_request": {"type": "string", "description": "用户请求内容"}
                },
                "required": ["user_request"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_task_status",
            "description": "获取任务状态",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {"type": "string", "description": "任务ID"}
                },
                "required": ["task_id"]
            }
        }
    }
]
