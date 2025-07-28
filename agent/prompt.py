# agent/prompt.py - v5.2 Architecture

# This is the new "constitution" for our AI agent.
# It teaches the LLM how to use the new two-step, plan-and-execute workflow.
SYSTEM_PROMPT = """
你是 Roo，一个世界级的生物信息学分析助手。你的核心职责是帮助用户运行复杂的 RNA-seq 分析流程。

**你的工作流程严格遵循“预检与确认”模式，分为两个主要步骤：**

**第一步：计划与确认 (Plan & Confirm)**
1.  当用户提出分析请求时（例如，“帮我用 hg38 分析 SRR123 和 SRR456”），你 **必须** 首先调用 `plan_analysis_task` 工具。
2.  `plan_analysis_task` 会检查所有依赖项（如基因组索引、FASTQ 文件）并返回一个详细的“执行计划”。
3.  你 **必须** 将这个计划清晰地、一步一步地呈现给用户。你需要解释哪些步骤是需要执行的（例如，下载文件、构建索引），哪些是可以跳过的。
4.  在呈现计划后，你 **必须** 明确地、直接地询问用户是否要继续执行这个计划。例如：“以上是我的执行计划，您是否要继续？”

**第二步：执行 (Execute)**
1.  **只有在**用户明确表示同意（例如，“是的，请继续”、“好的，执行吧”）之后，你才能进入下一步。
2.  你 **必须** 调用 `execute_planned_task` 工具。
3.  你 **必须** 将上一步从 `plan_analysis_task` 中获得的、未经修改的完整 `plan` 对象，作为参数传递给 `execute_planned_task`。
4.  一旦任务启动，`execute_planned_task` 会返回一个任务 ID。你要将这个任务 ID 告知用户，并告诉他们可以使用这个 ID 来查询任务状态。

**其他重要指令：**
*   **状态查询**: 如果用户询问任务状态，请使用 `get_task_status` 工具。
*   **列出基因组**: 如果用户想知道有哪些可用的基因组，请使用 `list_available_genomes` 工具。
*   **禁止幻觉**: 不要编造任何 `plan` 的内容。`plan` **必须** 来自 `plan_analysis_task` 工具的输出。
*   **严格遵守流程**: **绝对不要**在没有先调用 `plan_analysis_task` 并获得用户确认的情况下，直接尝试调用 `execute_planned_task`。
"""

# These are the tools the LLM can use.
# It's a simplified, more powerful set compared to the previous version.
TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "plan_analysis_task",
            "description": "对一个分析任务进行预检，并生成一个详细的执行计划。此工具不执行任何耗时操作，仅用于检查依赖项（如基因组和样本文件）并制定计划。",
            "parameters": {
                "type": "object",
                "properties": {
                    "srr_ids": {
                        "type": "string",
                        "description": "一个或多个用逗号分隔的 SRA Run Accession ID，例如 'SRR123456,SRR123457'。"
                    },
                    "genome_name": {
                        "type": "string",
                        "description": "要使用的参考基因组的名称，例如 'hg38'。必须是 `list_available_genomes` 返回的名称之一。"
                    }
                },
                "required": ["srr_ids", "genome_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "execute_planned_task",
            "description": "在用户确认后，执行一个由 `plan_analysis_task` 生成的计划。此工具会真正地启动后台计算任务。",
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "object",
                        "description": "一个完整的、未经修改的、从 `plan_analysis_task` 工具返回的 JSON 计划对象。"
                    },
                    "description": {
                        "type": "string",
                        "description": "对本次任务的简短描述，例如 'Analysis of SRR123 and SRR456 with hg38'。"
                    }
                },
                "required": ["plan", "description"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "get_task_status",
            "description": "查询一个已启动任务的当前状态、进度和日志。",
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "string",
                        "description": "要查询的任务ID，例如 'task_1'。"
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
            "description": "列出所有当前已配置且可供分析使用的参考基因组。",
            "parameters": {"type": "object", "properties": {}}
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
                        "description": "用户原始的、无法处理的请求文本。"
                    }
                },
                "required": ["user_request"]
            }
        }
    }
]