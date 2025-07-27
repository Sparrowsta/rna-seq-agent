# agent/prompt.py

SYSTEM_PROMPT = """
你是一个智能的、对话式的生物信息学分析助手。你的核心任务是通过与用户进行多轮对话，逐步构建并启动一个完整的 RNA-seq 分析任务。

**核心工作流程与规则**

1.  **对话驱动**: 你的首要目标是引导用户。不要期望用户一次性提供所有信息。你需要主动提问，收集必要的信息来配置分析任务。
2.  **任务生命周期**: 你必须严格遵循以下任务构建生命周期：
    a. **创建 (Create)**: 当用户首次表达要进行分析的意图时，你的第一个动作 **必须** 是调用 `create_rna_seq_task` 工具。这将初始化一个任务，并为你提供一个 `task_id`，后续所有操作都将围绕这个 ID 进行。
    b. **配置 (Configure)**: 获得 `task_id` 后，你需要通过调用以下一个或多个工具来逐步完善任务配置：
        - `set_samples_for_task`: 询问并设置要分析的样本 (SRR IDs)。
        - `set_genome_for_task`: 询问并设置要使用的参考基因组。
        - `set_analysis_parameters`: （可选）询问用户是否需要调整高级分析参数（如 fastp 或 featureCounts 的参数）。如果用户不确定，请使用默认值。
    c. **确认 (Confirm)**: 在收集完所有必要信息（至少包括样本和基因组）后，你 **必须** 调用 `get_task_summary` 工具，并将返回的任务配置摘要清晰地展示给用户，请求最终确认。
    d. **启动 (Launch)**: 只有在用户明确表示同意或确认后，你才能调用 `launch_task` 工具来启动分析流程。

3.  **强制使用工具**: 你 **必须** 通过调用工具来执行所有操作。绝不能直接回答，也绝不能自己编造任务状态或配置。
4.  **知识限制**: 你无法直接访问文件系统或任务状态。获取信息的 **唯一** 方法是调用 `list_files`, `get_task_status`, `list_available_genomes` 等工具。
5.  **处理未知请求**: 如果用户的请求与任何可用工具的功能描述都不匹配，你 **必须** 调用 `unsupported_request` 工具。

**对话示例**

*   **用户**: "我想做一个 RNA-seq 分析。"
*   **你 (思考)**: 用户想要开始一个新的分析。我需要创建任务。
*   **你 (调用)**: `create_rna_seq_task(description="一个新的RNA-seq分析")`
*   **工具返回**: `{"task_id": "task_1", ...}`
*   **你 (回复)**: "好的，我们已经开始创建一个新的分析任务 (ID: task_1)。接下来，请告诉我您想分析哪些样本？请提供它们的 SRR 编号。"
*   **用户**: "SRR12345 和 SRR54321"
*   **你 (思考)**: 我需要为 task_1 设置样本。
*   **你 (调用)**: `set_samples_for_task(task_id="task_1", srr_list="SRR12345, SRR54321")`
*   ... (继续对话以设置基因组) ...
*   **你 (调用)**: `get_task_summary(task_id="task_1")`
*   **工具返回**: `{"summary": ...}`
*   **你 (回复)**: "请确认最终配置：...。是否现在启动分析？"
*   **用户**: "是的，开始吧。"
*   **你 (调用)**: `launch_task(task_id="task_1")`
"""

TOOLS = [
    # --- 1. 任务构建核心流程 ---
    {
        "type": "function",
        "function": {
            "name": "create_rna_seq_task",
            "description": "当用户首次表达想要运行一个新的分析流程时，必须首先调用此工具。它会初始化一个任务并返回一个唯一的 task_id。",
            "parameters": {
                "type": "object",
                "properties": {
                    "description": {
                        "type": "string",
                        "description": "根据用户的请求，为任务生成一个简短的描述。例如：'分析人类细胞系的RNA-seq数据'。"
                    }
                },
                "required": ["description"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "set_samples_for_task",
            "description": "为一个已经创建的任务设置要分析的样本。你需要提供任务ID和样本的SRR列表。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'create_rna_seq_task' 返回的任务唯一标识符。"
                    },
                    "srr_list": {
                        "type": "string",
                        "description": "一个包含一个或多个 SRR (Sequence Read Archive) 运行编号的字符串，可以由逗号或空格分隔。例如: 'SRR12345, SRR67890'"
                    }
                },
                "required": ["task_id", "srr_list"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "set_genome_for_task",
            "description": "为一个已经创建的任务设置参考基因组。你需要提供任务ID和基因组名称。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'create_rna_seq_task' 返回的任务唯一标识符。"
                    },
                    "genome_name": {
                        "type": "string",
                        "description": "要用于分析的基因组的名称。可以使用 'list_available_genomes' 工具查看可用选项。"
                    }
                },
                "required": ["task_id", "genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "set_analysis_parameters",
            "description": "为一个已经创建的任务设置特定的分析工具参数。这是一个高级选项，仅当用户明确要求修改时使用。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'create_rna_seq_task' 返回的任务唯一标识符。"
                    },
                    "tool_name": {
                        "type": "string",
                        "enum": ["fastp", "featureCounts"],
                        "description": "要配置参数的工具名称。"
                    },
                    "params": {
                        "type": "object",
                        "description": "一个包含参数名和值的字典。例如: '{\"fastp_q_val\": 25, \"fc_is_paired_end\": false}'"
                    }
                },
                "required": ["task_id", "tool_name", "params"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_task_summary",
            "description": "在收集完所有必要信息后，调用此工具来获取任务的最终配置摘要，以便向用户进行确认。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'create_rna_seq_task' 返回的任务唯一标识符。"
                    }
                },
                "required": ["task_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "launch_task",
            "description": "在用户最终确认了任务配置后，调用此工具来启动实际的 Nextflow 分析流程。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'create_rna_seq_task' 返回并经过用户确认的任务唯一标识符。"
                    }
                },
                "required": ["task_id"]
            }
        }
    },
    # --- 2. 状态与文件管理 ---
    {
        "type": "function",
        "function": {
            "name": "get_task_status",
            "description": "当用户想要查询一个已经提交或正在运行的任务的状态时，调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "任务的唯一标识符，例如 'task_1' 或 'download_hg38'。"
                    }
                },
                "required": ["task_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": "列出指定路径下的文件和目录。路径是相对于 'data' 目录的。用于检查流程输出或已下载的文件。",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "要查看的目录路径，相对于 'data' 目录。例如，要查看 'data/results'，请提供 'results'。如果省略，则默认为 'data' 的根目录。"
                    }
                },
                "required": []
            }
        }
    },
    # --- 3. 基因组管理 ---
    {
        "type": "function",
        "function": {
            "name": "list_available_genomes",
            "description": "当用户想要查询当前服务器上有哪些可用的基因组时调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {},
                "required": []
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "add_genome_to_config",
            "description": "将一个新的基因组条目添加到配置文件中。这个操作只更新配置，不执行下载。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {"type": "string", "description": "为该基因组指定一个简短且唯一的名称，例如 'hg19'。"},
                    "species": {"type": "string", "description": "该基因组所属的物种，例如 'human'。"},
                    "version": {"type": "string", "description": "该基因组的版本号，例如 'hg19'。"},
                    "fasta_url": {"type": "string", "description": "指向 FASTA 文件的完整 URL。"},
                    "gtf_url": {"type": "string", "description": "指向 GTF 文件的完整 URL。"}
                },
                "required": ["genome_name", "species", "version", "fasta_url", "gtf_url"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "download_genome_files",
            "description": "为一个已经存在于配置文件中的基因组启动后台下载任务。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {
                        "type": "string",
                        "description": "要下载的基因组的名称，这个名称必须已经存在于配置文件中。"
                    }
                },
                "required": ["genome_name"]
            }
        }
    },
    # --- 4. 兜底工具 ---
    {
        "type": "function",
        "function": {
            "name": "unsupported_request",
            "description": "当用户的请求与任何其他可用工具的功能都不匹配时，必须调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {
                    "user_request": {
                        "type": "string",
                        "description": "用户的原始请求字符串。"
                    }
                },
                "required": ["user_request"]
            }
        }
    }
]