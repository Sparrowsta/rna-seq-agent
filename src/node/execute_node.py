import json
import subprocess
import asyncio
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any
from ..state import AgentState

async def execute_node(state: AgentState) -> Dict[str, Any]:
    """执行节点 - 构建和执行Nextflow命令"""
    print(f"\n{'='*60}")
    print(f"🚀 **RNA-seq分析执行**")
    print(f"{'='*60}")
    
    # 获取配置
    nextflow_config = state.nextflow_config or {}
    print(f"📊 **分析配置:**")
    for key, value in nextflow_config.items():
        print(f"   {key}: {value}")
    
    # 生成运行时配置文件
    print(f"\n📝 **生成运行时配置...**")
    config_result = await generate_runtime_config(nextflow_config)
    
    if not config_result["success"]:
        return {
            "nextflow_command": "",
            "execution_status": "failed",
            "execution_output": f"配置生成失败: {config_result['error']}",
            "execution_result": {"success": False, "error": config_result["error"]},
            "response": "分析执行失败：配置生成错误",
            "status": "failed"
        }
    
    # 构建Nextflow命令
    print(f"\n🔧 **构建Nextflow命令...**")
    nextflow_command = build_nextflow_command(nextflow_config)
    print(f"📋 命令: {nextflow_command}")
    
    # 执行Nextflow流水线
    print(f"\n⚡ **执行Nextflow流水线...**")
    print(f"🕐 开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    execution_result = await execute_nextflow_pipeline(nextflow_command)
    
    print(f"🕐 结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📊 **执行结果:** {'✅ 成功' if execution_result['success'] else '❌ 失败'}")
    
    # 生成响应消息
    if execution_result["success"]:
        response_msg = f"""🎉 **RNA-seq分析执行成功！**

📋 **执行摘要:**
   - 基因组版本: {nextflow_config.get('genome_version', 'unknown')}
   - 分析工具链: {nextflow_config.get('qc_tool', 'unknown')}-{nextflow_config.get('align_tool', 'unknown')}-{nextflow_config.get('quant_tool', 'unknown')}
   - 执行时长: {execution_result.get('duration', 'unknown')}
   - 输出目录: data/results/

💡 **下一步:** 查看 data/results/ 目录中的分析结果"""
    else:
        response_msg = f"""❌ **RNA-seq分析执行失败**

🔍 **错误信息:**
{execution_result.get('error', '未知错误')}

💡 **建议:** 检查配置参数和数据文件完整性"""
    
    return {
        "nextflow_command": nextflow_command,
        "execution_status": "completed" if execution_result["success"] else "failed",
        "execution_output": execution_result.get("output", ""),
        "execution_result": execution_result,
        "response": response_msg,
        "status": "execute"
    }

async def generate_runtime_config(nextflow_config: Dict[str, Any]) -> Dict[str, Any]:
    """生成运行时配置文件"""
    try:
        config_dir = Path("/config")
        config_dir.mkdir(exist_ok=True)
        
        # 创建运行时配置
        runtime_config = {
            "timestamp": datetime.now().isoformat(),
            "analysis_id": f"rna_seq_{int(time.time())}",
            "nextflow_params": nextflow_config
        }
        
        # 保存配置文件
        config_file = config_dir / "runtime_config.json"
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(runtime_config, f, indent=2, ensure_ascii=False)
        
        print(f"✅ 配置文件已保存: {config_file}")
        return {"success": True, "config_file": str(config_file)}
        
    except Exception as e:
        print(f"❌ 配置生成失败: {e}")
        return {"success": False, "error": str(e)}

def build_nextflow_command(nextflow_config: Dict[str, Any]) -> str:
    """构建Nextflow命令"""
    # 基础命令 - 从data目录执行根目录的main.nf
    cmd_parts = ["nextflow", "run", "/main.nf"]
    
    # 添加参数
    if nextflow_config.get("genome_version"):
        cmd_parts.extend(["--genome_version", nextflow_config["genome_version"]])
    
    if nextflow_config.get("qc_tool"):
        cmd_parts.extend(["--qc_tool", nextflow_config["qc_tool"]])
    
    if nextflow_config.get("align_tool"):
        cmd_parts.extend(["--align_tool", nextflow_config["align_tool"]])
    
    if nextflow_config.get("quant_tool"):
        cmd_parts.extend(["--quant_tool", nextflow_config["quant_tool"]])
    
    # 样本配对信息 - Agent分析的结果，包含完整的文件路径信息
    sample_groups = nextflow_config.get("sample_groups", [])
    if sample_groups:
        import json
        # 将样本配对信息转换为JSON字符串传递给Nextflow
        sample_groups_json = json.dumps(sample_groups, separators=(',', ':'))
        cmd_parts.extend(["--sample_groups", f"'{sample_groups_json}'"])
    
    # 基因组和索引管理
    if nextflow_config.get("run_download_genome"):
        cmd_parts.extend(["--run_download_genome", str(nextflow_config["run_download_genome"]).lower()])
    
    if nextflow_config.get("run_build_star_index"):
        cmd_parts.extend(["--run_build_star_index", str(nextflow_config["run_build_star_index"]).lower()])
    
    # 工作目录设置 - 使用相对路径
    cmd_parts.extend(["-work-dir", "work"])
    # 生成清理选项(可选)
    cmd_parts.append("-resume")  # 支持断点续传
    
    return " ".join(cmd_parts)

async def execute_nextflow_pipeline(command: str) -> Dict[str, Any]:
    """执行Nextflow流水线"""
    start_time = time.time()
    
    try:
        print(f"🔄 执行命令: {command}")
        
        # 在Docker容器环境中执行，当前工作目录应该已经是/data
        print(f"🚀 启动Nextflow执行...")
        process = await asyncio.create_subprocess_shell(
            command,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
            cwd="."
        )
        
        # 实时读取输出
        output_lines = []
        while True:
            line = await process.stdout.readline()
            if not line:
                break
            
            line_text = line.decode('utf-8').strip()
            output_lines.append(line_text)
            
            # 实时显示关键信息
            if any(keyword in line_text.lower() for keyword in 
                   ['completed', 'failed', 'error', 'submitted', 'cached']):
                print(f"   📋 {line_text}")
        
        # 等待进程完成
        await process.wait()
        
        # 计算执行时间
        duration = time.time() - start_time
        duration_str = f"{duration:.1f}秒"
        
        full_output = "\n".join(output_lines)
        
        if process.returncode == 0:
            print(f"✅ Nextflow执行成功 (耗时: {duration_str})")
            return {
                "success": True,
                "output": full_output,
                "duration": duration_str,
                "return_code": process.returncode,
                "mode": "real"
            }
        else:
            print(f"❌ Nextflow执行失败 (返回码: {process.returncode})")
            return {
                "success": False,
                "output": full_output,
                "duration": duration_str,
                "return_code": process.returncode,
                "error": f"Nextflow进程失败，返回码: {process.returncode}",
                "mode": "real"
            }
    
    except Exception as e:
        duration = time.time() - start_time
        error_msg = f"执行异常: {str(e)}"
        print(f"❌ {error_msg}")
        
        return {
            "success": False,
            "output": "",
            "duration": f"{duration:.1f}秒",
            "error": error_msg,
            "mode": "error"
        }