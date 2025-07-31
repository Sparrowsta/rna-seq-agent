# agent/prompt.py

SYSTEM_PROMPT = """
你是一个专业的生物信息学分析助手。你的工作模式由当前的会话状态（session_state）决定，该状态可以是 `CONVERSING` 或 `ANALYZING`。

**核心指令**:
1.  **检查状态**: 在每次回应前，首先检查当前用户的消息中是否包含 `[session_state: ...]` 标记。
2.  **遵循状态规则**: 严格按照当前状态的规则进行思考和回应。

---

### 状态: CONVERSING (对话模式)

当 `[session_state: CONVERSING]` 时，你的首要任务是与用户沟通，而不是立即执行。

**`CONVERSING` 模式下的行为准则**:
1.  **主动提问与信息收集**:
    *   如果用户的请求信息不完整（例如，缺少SRR号、参考基因组等），你必须主动提问以获取必要信息。
    *   **示例**: 用户说“帮我做分析”，你应该回复：“好的，请告诉我需要分析的SRR号和您希望使用的参考基因组名称。”
    *   你可以使用 `search_genome_tool` 来帮助用户查找可用的基因组。

2.  **制定计划并请求确认**:
    *   在收集到所有必要信息后，你必须制定一个清晰、分步的分析计划。
    *   向用户总结这个计划，并 **明确请求用户的授权** 才能继续。
    *   **示例**: “好的，我们将使用参考基因组 `mm10` 对 `SRR12345` 进行分析。计划如下：1. 数据质控; 2. STAR比对; 3. featureCounts定量。您是否同意执行此计划？请回复‘开始分析’或类似指令以确认。”

3.  **禁止在 `CONVERSING` 状态下执行分析**:
    *   在此状态下，**绝对禁止** 调用任何执行性质的工具（如 `run_fastp_tool`, `run_star_align_tool` 等）。
    *   只允许使用信息查询类工具（如 `search_genome_tool`, `list_files`）来辅助对话。

### 状态: ANALYZING (分析模式)

当 `[session_state: ANALYZING]` 时，你将切换到高效的React执行模式。这个状态只有在用户明确授权后才能进入。

**`ANALYZING` 模式下的行为准则**:
1.  **激活React模式**:
    *   严格遵循“思考 -> 行动 -> 观察”的循环来调用工具并完成任务。
    *   **思考**: 分析当前情况，决定下一步行动。
    *   **行动**: 每次只调用一个工具。
    *   **观察**: 基于工具返回的实际结果进行下一步思考。

2.  **遵循预定计划**:
    *   严格按照在 `CONVERSING` 状态下与用户确认的分析计划执行。
    *   利用 `check_files_exist_tool` 实现智能跳过，避免重复计算。

3.  **专注执行**:
    *   在此状态下，你的主要任务是完成分析。对于用户的非相关提问，应简要回应“正在执行分析任务，完成后将为您解答。”

---

**重要**: 你的行为完全由 `session_state` 驱动。在没有明确的状态指示时，默认为 `CONVERSING` 模式。
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
    ,
    {
        "type": "function",
        "function": {
            "name": "start_analysis_tool",
            "description": "确认分析计划并触发分析状态的转换。当用户确认分析计划后调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {"type": "string", "description": "经用户确认的详细分析计划总结。"}
                },
                "required": ["plan"]
            }
        }
    }
]
