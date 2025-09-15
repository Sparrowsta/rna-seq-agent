import json
from typing import Dict, Any
from ..state import AgentState, PrepareResponse
from ..core import get_shared_llm
from ..prompts import PREPARE_NODE_PROMPT
from langgraph.prebuilt import create_react_agent
from ..tools import (
    scan_fastq_files,
    scan_system_resources,
    scan_genome_files,
    check_tool_availability,
    get_project_overview,
    write_params_file,
)
from ..logging_bootstrap import get_logger, log_llm_preview

logger = get_logger("rna.nodes.prepare")

def create_prepare_agent(detection_context: str = ""):
    """创建Prepare节点的智能配置Agent
    
    Args:
        detection_context: 检测数据和用户需求的上下文信息
    """
    llm = get_shared_llm()
    
    # 构建完整的prompt，结合基础prompt和上下文
    if detection_context:
        full_prompt = PREPARE_NODE_PROMPT + "\n\n## 当前分析数据\n" + detection_context
    else:
        full_prompt = PREPARE_NODE_PROMPT
    
    # 使用create_react_agent并挂载检测相关工具，允许按需再检测
    agent = create_react_agent(
        model=llm,
        tools=[
            scan_fastq_files,
            scan_system_resources,
            scan_genome_files,
            check_tool_availability,
            get_project_overview,
        ],
        prompt=full_prompt,  # 使用动态构建的完整prompt
        response_format=PrepareResponse
    )
    return agent


async def prepare_node(state: AgentState) -> Dict[str, Any]:
    """准备节点 - 专注于初始配置生成，基于用户需求和检测数据"""
    logger.info("开始智能配置分析...")
    
    # 获取核心信息
    detection_results = state.query_results or {}
    initial_requirements = state.user_requirements or {}
    
    if not detection_results:
        return {
            "success": False,
            "nextflow_config": {},
            "resource_config": {},
            "config_reasoning": "未获取到检测数据，无法进行智能配置分析",
            "response": "⚠️ 缺少检测数据，无法进行智能配置分析",
            "status": "failed"
        }
    
    # 使用LLM综合分析用户需求和检测数据
    try:
        logger.info("LLM正在基于用户需求分析检测数据并生成初始配置...")
        
        # 构建上下文（仅包含初始需求，不处理修改）
        context_parts = []
        if initial_requirements:
            context_parts.append(f"用户需求: {initial_requirements}")
        
        # 添加检测数据
        context_parts.append(json.dumps(detection_results, indent=2, ensure_ascii=False))
        
        detection_context = "\n".join(context_parts)
        
        # 将上下文传递给create_prepare_agent
        agent_executor = create_prepare_agent(detection_context)
        
        # 构建用户消息
        user_message = "请基于用户需求和检测数据生成最优初始配置"
        messages_input = {"messages": [{"role": "user", "content": user_message}]}
        
        result = await agent_executor.ainvoke(messages_input)
        structured_response = result.get("structured_response")
        # DEBUG: 记录结构化输出预览
        try:
            if structured_response:
                log_llm_preview(logger, "prepare", structured_response)
            else:
                # 结构化缺失时，记录关键字段预览
                log_llm_preview(logger, "prepare.raw", {"keys": list(result.keys())[:10]})
        except Exception:
            pass

        # 检查LLM响应并提取结果
        if structured_response:
            reasoning = structured_response.config_reasoning or "基于用户需求和检测数据的智能分析"

            nextflow_cfg = structured_response.nextflow_config or {}
            resource_params = structured_response.resource_config or {}

            logger.info("初始配置生成完成")

            # 参数版本化: 合并配置并写入版本化文件 (M1实现)
            try:
                # 合并 nextflow_config 和 resource_config
                merged_params = {
                    "nextflow_config": nextflow_cfg,
                    "resource_config": resource_params
                }
                
                # 写入参数版本化文件 (仅写入，不参与运行)
                param_file_path = write_params_file("prepare", merged_params, state)
                logger.info(f"Prepare参数版本化文件已写入: {param_file_path}")
                
            except Exception as e:
                logger.warning(f"参数版本化写入失败 (不影响配置生成): {e}")

            # 构建用户需求满足说明
            user_satisfaction_note = ""
            if initial_requirements:
                user_satisfaction_note = f"\n\n🎯 用户需求处理情况：\n📋 {initial_requirements}"

            return {
                "success": True,
                "nextflow_config": nextflow_cfg,
                "resource_config": resource_params,
                "config_reasoning": reasoning,
                "response": f"智能配置分析完成{user_satisfaction_note}\n\n💡 {reasoning}",
                "status": "success"
            }
        else:
            raise Exception("Agent未返回预期的结构化响应")
        
    except Exception as e:
        logger.error(f"配置生成失败: {str(e)}", exc_info=True)
        return {
            "success": False,
            "nextflow_config": {},
            "resource_config": {},
            "config_reasoning": f"配置生成失败: {str(e)}",
            "response": f"❌ 配置生成失败: {str(e)}",
            "status": "failed"
        }
