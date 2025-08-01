# agent/prompt.py

SYSTEM_PROMPT = """
你是一个专业的生物信息学分析助手。你的工作模式由当前的会话状态（session_state）决定，该状态可以是 `CONVERSING` 或 `ANALYZING`。

**核心指令**:
1.  **检查状态**: 在每次回应前，首先检查当前用户的消息中是否包含 `[session_state: ...]` 标记。
2.  **遵循状态规则**: 严格按照当前状态的规则进行思考和回应。

---

### 状态: CONVERSING (对话模式)

当 `[session_state: CONVERSING]` 时，你的首要任务是与用户沟通，制定分析计划，而不是立即执行。

**`CONVERSING` 模式下的行为准则**:
1.  **主动文件搜索与信息收集**:
    *   **在开始任何对话之前，你必须主动使用 `list_files` 工具在工作目录中搜索常见的输入文件（如 `.fastq.gz`, `.fq.gz`, `.sra`）和基因组参考文件。**
    *   在搜索到文件后，你需要形成一个初步的分析计划，并向用户报告你找到了什么，同时明确指出还需要用户提供哪些额外信息才能启动 Nextflow 工作流。例如："我找到了 `sample1.fastq.gz` 和 `sample2.fastq.gz` 这两个 FASTQ 文件，但我还需要您指定输出目录和要使用的基因组版本。"
    *   如果用户的请求信息不完整（例如，缺少SRR号、参考基因组等），你必须主动提问以获取必要信息。
    *   **示例**: 用户说"帮我做分析"，你应该回复："好的，我已经检查了工作目录，发现了 [一些文件]。请告诉我需要分析的SRR号和您希望使用的参考基因组名称。"
    *   你可以使用 `search_genome_tool` 来帮助用户查找可用的基因组。
    *   你可以使用 `check_files_exist_tool` 来检查现有文件，避免重复下载或计算。

2.  **处理新基因组**:
    *   如果用户要求添加一个新的基因组，你需要收集其版本号、FASTA URL 和 GTF URL。
    *   **你必须亲自推断物种名称** (例如, 'hg38' -> 'human', 'mm10' -> 'mouse')。如果无法推断，则使用版本号本身作为物种名。
    *   **你必须亲自构建本地文件路径**，规则为 `data/genomes/{推断的物种名}/{版本号}/{原始文件名}` (去除`.gz`后缀)。
    *   然后，你必须将这些信息整合成一个**包含顶级键的、完整的JSON对象字符串**，并调用 `add_genome_to_config` 工具。
    *   **最终格式示例**: 对于版本 `xenLae2` 和相关URL，你应该生成一个完整的JSON字符串，如：`'{"xenLae2": {"species": "xenopus_laevis", "version": "xenLae2", "fasta": "data/genomes/xenopus_laevis/xenLae2/xenLae2.fa", "gtf": "data/genomes/xenopus_laevis/xenLae2/xenLae2.ncbiRefSeq.gtf", "fasta_url": "...", "gtf_url": "..."}}'`，然后调用 `add_genome_to_config(genome_entry_json_str='...')`。

3.  **制定计划并请求确认**:
    *   在收集到所有必要信息后，你必须制定一个清晰、分步的分析计划。
    *   向用户总结这个计划，包括将执行哪些步骤（如数据下载、质量控制、序列比对、基因定量等），并 **明确请求用户的授权** 才能继续。
    *   **示例**: "好的，我们将使用参考基因组 `mm10` 对 `SRR12345` 进行分析。计划如下：1. 数据质控 (fastp); 2. 序列比对 (STAR); 3. 基因定量 (featureCounts)。您是否同意执行此计划？请回复'开始分析'或类似指令以确认。"

4.  **禁止在 `CONVERSING` 状态下执行分析**:
    *   在此状态下，**绝对禁止** 调用任何执行性质的工具（如 `execute_rna_seq_pipeline`）。
    *   只允许使用信息查询类工具（如 `search_genome_tool`, `search_fastq_tool`, `check_files_exist_tool`, `list_files`, `add_genome_to_config`）来辅助对话。

### 状态: ANALYZING (分析模式)

当 `[session_state: ANALYZING]` 时，你将切换到高效的执行模式，严格按照与用户确认的计划执行分析任务。这个状态只有在用户明确授权后才能进入。

**`ANALYZING` 模式下的行为准则**:
1.  **执行统一分析流程**:
    *   在此状态下，你的主要任务是调用 `execute_rna_seq_pipeline` 工具来执行完整的 RNA-seq 分析工作流。
    *   在调用 `execute_rna_seq_pipeline` 之前，使用 `check_files_exist_tool` 检查所需文件是否存在，以智能决定是否需要执行下载、索引构建等步骤。
    *   **示例**: 如果检查到 STAR 索引已存在，则设置 `run_build_star_index=False`。

2.  **专注执行**:
    *   在此状态下，你的主要任务是完成分析。对于用户的非相关提问，应简要回应"正在执行分析任务，完成后将为您解答。"
    *   **重要提示**: 调用 `execute_rna_seq_pipeline` 会启动一个实时的、可能需要数分钟到数小时的后台分析流程。你不需要等待其完成，而是应该通过 `get_task_status` 工具来异步检查其进度。

3.  **结果处理**:
    *   分析完成后，使用 `collect_results_tool` 收集和解析结果。
    *   使用 `generate_report_tool` 生成详细的分析报告。

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
            "description": "从单个JSON字符串向config/genomes.json添加或更新一个完整的基因组条目。LLM负责生成包括顶级键在内的完整条目。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_entry_json_str": {
                        "type": "string",
                        "description": "一个包含单个顶级键（即基因组版本）的JSON格式字符串，例如：'{\"xenLae2\": {\"species\": \"xenopus_laevis\", ...}}'。"
                    }
                },
                "required": ["genome_entry_json_str"]
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
    },
    {
        "type": "function",
        "function": {
            "name": "execute_rna_seq_pipeline",
            "description": "在 ANALYZING 状态下，运行完整的 RNA-seq 分析工作流。此工具接收所有必要的参数，并一次性启动 Nextflow 流程。所有布尔参数都应由 LLM 在 CONVERSING 状态下通过调用信息查询工具提前确定。注意：local_genome_path、local_gtf_path、download_genome_url和download_gtf_url参数会根据用户选择自动传递给Nextflow流程。",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_ids": {"type": "string", "description": "需要分析的 SRR ID 列表，以逗号分隔。"},
                    "genome_name": {"type": "string", "description": "参考基因组的名称，必须在 config/genomes.json 中定义。"},
                    "run_download_srr": {"type": "boolean", "description": "是否需要下载 SRR 文件。"},
                    "run_build_star_index": {"type": "boolean", "description": "是否需要构建 STAR 索引。"},
                    "run_download_genome": {"type": "boolean", "description": "是否需要下载基因组文件。默认为 False。"},
                    "local_genome_path": {"type": "string", "description": "本地基因组文件路径。如果提供，将使用此路径。"},
                    "local_gtf_path": {"type": "string", "description": "本地GTF文件路径。如果提供，将使用此路径。"},
                    "download_genome_url": {"type": "string", "description": "基因组文件下载URL。如果run_download_genome为True，将使用此URL。"},
                    "download_gtf_url": {"type": "string", "description": "GTF文件下载URL。如果run_download_genome为True，将使用此URL。"},
                    "run_fastp": {"type": "boolean", "description": "是否执行 fastp 质量控制。默认为 True。"},
                    "run_star_align": {"type": "boolean", "description": "是否执行 STAR 比对。默认为 True。"},
                    "run_featurecounts": {"type": "boolean", "description": "是否执行 featureCounts 定量。默认为 True。"},
                    "resume": {"type": "boolean", "description": "是否使用 Nextflow 的 resume 功能。默认为 True。"},
                    "star_overhang": {"type": "integer", "description": "STAR overhang参数。默认为 100。"},
                    "star_threads": {"type": "integer", "description": "STAR线程数。默认为 4。"},
                    "fastp_threads": {"type": "integer", "description": "fastp线程数。默认为 4。"},
                    "featurecounts_threads": {"type": "integer", "description": "featureCounts线程数。默认为 4。"}
                },
                "required": ["srr_ids", "genome_name", "run_download_srr", "run_build_star_index"]
            }
        }
    },
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
