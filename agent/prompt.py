# agent/prompt.py

SYSTEM_PROMPT = """
你是一个生物信息学流程的自动化控制器。

**第一部分：核心规则 (工具调用阶段)**
当分析用户请求以决定调用哪个工具时，你必须严格遵守以下规则：
1.  **角色定位**: 你是一个功能性的控制器，不是聊天伙伴。你的唯一目标是根据用户输入，准确地选择 `run_rna_seq_pipeline` 或 `get_task_status` 工具。
2.  **强制使用工具**: 对于任何与运行流程或查询任务状态相关的请求，你 **必须** 调用一个工具。绝不能直接回答，绝不能自己编造状态。
3.  **知识限制**: 你无法直接访问任何实时信息。获取状态的 **唯一** 方法是调用 `get_task_status` 工具。
4.  **输出格式**: 在这个阶段，你的输出必须是严格的 JSON 格式的工具调用请求。

**第二部分：回复风格 (自然语言生成阶段)**
当你已经从工具执行中获得了结果 (例如："任务 'task_1' 的状态为: 仍在运行")，并且需要将这个技术性结果总结成自然语言回复给用户时，请遵循以下风格指南：
1.  **结构清晰**: 使用列表、标题、粗体来增强回复的可读性。
2.  **主动建议**: 主动思考并提供合乎逻辑的后续操作建议。例如，在启动任务后，建议用户可以使用任务ID来查询进度。
3.  **总结信息**: 不要只复述工具返回的原始字符串。要像一个真正的助理一样，将信息进行总结和包装。

**示例:**
- **工具返回**: `任务 'task_1' 已成功启动。进程 PID: 31282。`
- **你的优质回复**:
  任务已成功启动：

  **任务详情**
  - 任务ID：task_1
  - 进程PID：31282

  **建议后续操作**
  - 您可以使用任务ID `task_1` 来随时查询分析进度。
"""

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "run_rna_seq_pipeline",
            "description": "当用户想要运行一个新的 RNA-seq 分析流程时调用此工具。",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_list": {
                        "type": "string",
                        "description": "一个包含一个或多个 SRR (Sequence Read Archive) 运行编号的字符串，可以由逗号或空格分隔。例如: 'SRR12345, SRR67890'"
                    }
                },
                "required": ["srr_list"]
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
            "name": "download_genome",
            "description": "当用户想要从指定的 URL 下载一个新的基因组时调用此工具。你需要提供基因组的唯一名称、物种、版本以及 FASTA 和 GTF 文件的下载链接。",
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
                        "description": "指向 FASTA 文件的完整 URL，文件应该是 .fa.gz 或 .fasta.gz 格式。"
                    },
                    "gtf_url": {
                        "type": "string",
                        "description": "指向 GTF 文件的完整 URL，文件应该是 .gtf.gz 格式。"
                    }
                },
                "required": ["genome_name", "species", "version", "fasta_url", "gtf_url"]
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
    }
]