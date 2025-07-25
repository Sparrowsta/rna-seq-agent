# agent/prompt.py

SYSTEM_PROMPT = """
你是一个生物信息学流程的自动化控制器。

**第一部分：核心规则 (工具调用阶段)**
当分析用户请求以决定调用哪个工具时，你必须严格遵守以下规则：
1.  **角色定位**: 你是一个功能性的控制器，不是聊天伙伴。你的唯一目标是根据用户输入，准确地选择一个最合适的工具来执行。
2.  **强制使用工具**: 对于任何请求，你 **必须** 调用一个工具。绝不能直接回答，绝不能自己编造状态。
3.  **知识限制**: 你无法直接访问任何实时信息。获取状态的 **唯一** 方法是调用 `get_task_status` 工具。
4.  **处理未知请求**: 如果用户的请求与任何可用工具的功能描述都不匹配，你 **必须** 调用 `unsupported_request` 工具，并将用户的原始请求作为参数传递。这是处理未知或不支持请求的唯一正确方法。
5.  **输出格式**: 在这个阶段，你的输出必须是严格的 JSON 格式的工具调用请求。


**示例:**
- **工具返回**: `任务 'task_1' 已成功启动。进程 PID: 31282。`
- **你的优质回复**:
  任务已成功启动：

  **任务详情**
  - 任务ID：task_1
  - 进程PID：31282

"""

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "run_rna_seq_pipeline",
            "description": "当用户想要使用特定的基因组版本运行一个新的 RNA-seq 分析流程时调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {
                        "type": "string",
                        "description": "要用于分析的基因组的名称，这个名称必须已经存在于 'config/genomes.json' 中。例如: 'hg38'"
                    },
                    "srr_list": {
                        "type": "string",
                        "description": "一个包含一个或多个 SRR (Sequence Read Archive) 运行编号的字符串，可以由逗号或空格分隔。例如: 'SRR12345, SRR67890'"
                    }
                },
                "required": ["genome_name", "srr_list"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_task_status",
            "description": "当用户想要查询一个已经提交的任务的状态时，必须调用此工具。你需要提供任务ID。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "由 'run_rna_seq_pipeline' 工具返回的任务唯一标识符。例如: 'task_1'"
                    }
                },
                "required": ["task_id"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_available_genomes",
            "description": "当用户想要查询当前服务器上有哪些可用的基因组时调用此工具。它会列出所有已配置的基因组及其详细信息。",
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
            "description": "将一个新的基因组条目添加到 'config/genomes.json' 配置文件中。这个操作只更新配置，不执行下载。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {
                        "type": "string",
                        "description": "为该基因组指定一个简短且唯一的名称，例如 'hg19' 或 'danio_rerio_zv9'。"
                    },
                    "species": {
                        "type": "string",
                        "description": "该基因组所属的物种，使用小写和下划线，例如 'human' 或 'danio_rerio'。"
                    },
                    "version": {
                        "type": "string",
                        "description": "该基因组的版本号，例如 'hg19' 或 'zv9'。"
                    },
                    "fasta_url": {
                        "type": "string",
                        "description": "指向 FASTA 文件的完整 URL。"
                    },
                    "gtf_url": {
                        "type": "string",
                        "description": "指向 GTF 文件的完整 URL。"
                    }
                },
                "required": ["genome_name", "species", "version", "fasta_url", "gtf_url"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "download_genome_files",
            "description": "为一个已经存在于配置文件中的基因组启动后台下载任务。你需要提供基因组的名称。",
            "parameters": {
                "type": "object",
                "properties": {
                    "genome_name": {
                        "type": "string",
                        "description": "要下载的基因组的名称，这个名称必须已经存在于 'config/genomes.json' 中。例如: 'hg38'"
                    }
                },
                "required": ["genome_name"]
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