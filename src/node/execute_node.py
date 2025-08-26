import json
import subprocess
import asyncio
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional
from ..state import AgentState

async def generate_nextflow_config(resource_config: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """生成动态的nextflow.config文件，包含资源配置"""
    try:
        config_dir = Path("/config")
        config_dir.mkdir(exist_ok=True)
        
        # 基础配置模板
        config_content = """// 动态生成的Nextflow配置文件
// 基于LLM智能资源分配

// 进程资源配置
process {
    // 默认配置
    cpus = 1
    memory = '2 GB'
    
    // 错误处理
    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries = 3
    maxErrors = '-1'
    
"""

        # 添加LLM智能分配的资源配置
        if resource_config:
            config_content += "    // LLM智能资源分配\n"
            for process_name, config in resource_config.items():
                cpus = config.get('cpus', 1)
                memory = config.get('memory', '2 GB')
                reasoning = config.get('reasoning', '默认配置')
                
                config_content += f"""    withName: '{process_name}' {{
        cpus = {cpus}
        memory = '{memory}'
        // {reasoning}
    }}
    
"""
        else:
            # 使用默认的硬编码配置
            config_content += """    // 默认资源配置（LLM未生成资源分配）
    withName: 'prepare_star_index' {
        cpus = 8
        memory = '32 GB'
        // 索引构建CPU密集
    }
    
    withName: 'run_quality_control' {
        cpus = 8
        memory = '16 GB'
        // 质控处理
    }
    
    withName: 'run_alignment' {
        cpus = 8
        memory = '32 GB'
        // 序列比对
    }
    
    withName: 'run_quantification' {
        cpus = 8
        memory = '16 GB'
        // 基因定量
    }
    
    withName: 'download_genome_fasta' {
        cpus = 2
        memory = '4 GB'
        // FASTA下载
    }
    
    withName: 'download_genome_gtf' {
        cpus = 2
        memory = '4 GB'
        // GTF下载
    }
    
"""

        # 添加执行器和报告配置
        config_content += """}

// 执行配置
executor {
    name = 'local'
}

// 报告配置
report {
    enabled = true
    file = 'results/nextflow/execution_report.html'
    overwrite = true
}

timeline {
    enabled = true
    file = 'results/nextflow/execution_timeline.html'
    overwrite = true
}

trace {
    enabled = true
    file = 'results/nextflow/execution_trace.txt'
    overwrite = true
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,read_bytes,write_bytes'
}
"""

        # 写入配置文件
        config_file = config_dir / "nextflow.config"
        with open(config_file, 'w', encoding='utf-8') as f:
            f.write(config_content)
        
        # 计算配置统计
        total_processes = len(resource_config) if resource_config else 6
        total_cpus = sum(config.get('cpus', 1) for config in resource_config.values()) if resource_config else 'default'
        
        print(f"✅ Nextflow配置文件已生成: {config_file}")
        print(f"📊 资源配置: {total_processes}个进程，总CPU分配: {total_cpus}")
        
        return {
            "success": True, 
            "config_file": str(config_file),
            "process_count": total_processes,
            "total_cpus": total_cpus
        }
        
    except Exception as e:
        print(f"❌ Nextflow配置生成失败: {e}")
        return {"success": False, "error": str(e)}

async def generate_runtime_config(nextflow_config: Dict[str, Any], resource_config: Optional[Dict[str, Dict[str, Any]]] = None) -> Dict[str, Any]:
    """生成运行时配置文件"""
    try:
        config_dir = Path("/config")
        config_dir.mkdir(exist_ok=True)
        
        # 明确处理resource_config的None情况
        if resource_config is None:
            resource_config = {}
        
        # 创建运行时配置
        runtime_config = {
            "timestamp": datetime.now().isoformat(),
            "analysis_id": f"rna_seq_{int(time.time())}",
            "nextflow_params": nextflow_config,
            "resource_config": resource_config
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
    
    # 明确指定配置文件路径
    cmd_parts.extend(["-c", "/config/nextflow.config"])
    
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
    
    # 明确传递下载和构建参数（无论true还是false）
    cmd_parts.extend(["--run_download_genome", str(nextflow_config.get("run_download_genome", False)).lower()])
    cmd_parts.extend(["--run_build_star_index", str(nextflow_config.get("run_build_star_index", False)).lower()])
    cmd_parts.extend(["--run_build_hisat2_index", str(nextflow_config.get("run_build_hisat2_index", False)).lower()])
    
    # 工作目录设置 - 使用相对路径
    cmd_parts.extend(["-work-dir", "work"])
    # 生成清理选项(可选)
    cmd_parts.append("-resume")  # 支持断点续传
    
    return " ".join(cmd_parts)

async def execute_nextflow_pipeline(command: str) -> Dict[str, Any]:
    """执行Nextflow流水线，包含SSL错误重试机制"""
    
    max_retries = 3
    ssl_error_keywords = [
        "ssl routines::unexpected eof while reading",
        "cannot download nextflow required file", 
        "make sure you can connect to the internet",
        "curl: (35) error:0a000126",
        "downloading nextflow dependencies"
    ]
    
    for attempt in range(max_retries):
        start_time = time.time()
        
        if attempt > 0:
            print(f"🔄 第 {attempt + 1} 次重试 (SSL错误重试)...")
            # 清理可能的残留进程和锁文件
            cleanup_cmd = "rm -rf .nextflow* work/.nextflow* || true"
            await asyncio.create_subprocess_shell(cleanup_cmd)
            await asyncio.sleep(2)  # 等待清理完成
        
        try:
            print(f"🔄 执行命令: {command}")
            print(f"🚀 启动Nextflow执行...")
            
            process = await asyncio.create_subprocess_shell(
                command,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
                cwd="."
            )
            
            # 实时读取输出
            output_lines = []
            ssl_error_detected = False
            
            # 类型断言确保stdout不为None
            assert process.stdout is not None, "stdout should not be None when PIPE is specified"
            
            while True:
                line = await process.stdout.readline()
                if not line:
                    break
                
                line_text = line.decode('utf-8').strip()
                output_lines.append(line_text)
                
                # 检查SSL错误
                if any(keyword in line_text.lower() for keyword in ssl_error_keywords):
                    print(f"⚠️ 检测到SSL网络错误: {line_text}")
                    ssl_error_detected = True
                    
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
            
            # 如果检测到SSL错误且还有重试机会
            if ssl_error_detected and attempt < max_retries - 1:
                print(f"❌ SSL错误导致执行失败，准备重试...")
                continue
            
            # 如果SSL错误但已到最后一次
            if ssl_error_detected and attempt == max_retries - 1:
                print(f"❌ SSL错误，已达最大重试次数 ({max_retries})")
                return {
                    "success": False,
                    "output": full_output,
                    "duration": duration_str,
                    "return_code": process.returncode,
                    "error": f"SSL网络错误，重试 {max_retries} 次后仍然失败",
                    "mode": "ssl_retry_failed"
                }
            
            # 正常执行结果判断
            if process.returncode == 0:
                success_msg = "✅ Nextflow执行成功"
                if attempt > 0:
                    success_msg += f" (第 {attempt + 1} 次重试成功)"
                success_msg += f" (耗时: {duration_str})"
                print(success_msg)
                
                return {
                    "success": True,
                    "output": full_output,
                    "duration": duration_str,
                    "return_code": process.returncode,
                    "mode": "success" if attempt == 0 else f"success_after_{attempt+1}_retries"
                }
            else:
                # 非SSL错误，不需要重试
                print(f"❌ Nextflow执行失败 (返回码: {process.returncode})")
                return {
                    "success": False,
                    "output": full_output,
                    "duration": duration_str,
                    "return_code": process.returncode,
                    "error": f"Nextflow进程失败，返回码: {process.returncode}",
                    "mode": "failed"
                }
        
        except Exception as e:
            duration = time.time() - start_time
            if attempt < max_retries - 1:
                print(f"❌ 执行异常，准备重试: {str(e)}")
                continue
            else:
                error_msg = f"执行异常: {str(e)}"
                print(f"❌ {error_msg}")
                return {
                    "success": False,
                    "output": "",
                    "duration": f"{duration:.1f}秒",
                    "error": error_msg,
                    "mode": "exception"
                }
    
    # 理论上不会到这里，但保险起见
    return {
        "success": False,
        "output": "所有重试尝试均失败",
        "duration": "unknown",
        "error": f"经过 {max_retries} 次重试仍然失败",
        "mode": "all_retries_failed"
    }

async def execute_node(state: AgentState) -> Dict[str, Any]:
    """执行节点 - 构建和执行Nextflow命令"""
    print(f"\n{'='*60}")
    print(f"🚀 **RNA-seq分析执行**")
    print(f"{'='*60}")
    
    # 获取配置
    nextflow_config = state.nextflow_config or {}
    resource_config = state.resource_config or {}
    
    print(f"📊 **分析配置:**")
    for key, value in nextflow_config.items():
        print(f"   {key}: {value}")
        
    if resource_config:
        print(f"🖥️ **资源配置:**")
        for process, config in resource_config.items():
            print(f"   {process}: {config.get('cpus')}核, {config.get('memory')}")
    
    # 生成动态的nextflow.config文件
    print(f"\n⚙️ **生成Nextflow配置文件...**")
    config_generation_result = await generate_nextflow_config(resource_config)
    
    if not config_generation_result["success"]:
        return {
            "nextflow_command": "",
            "execution_status": "failed",
            "execution_output": f"Nextflow配置生成失败: {config_generation_result['error']}",
            "execution_result": {"success": False, "error": config_generation_result["error"]},
            "response": "分析执行失败：Nextflow配置文件生成错误",
            "status": "failed"
        }
    
    # 生成运行时配置文件
    print(f"\n📝 **生成运行时配置...**")
    runtime_result = await generate_runtime_config(nextflow_config, resource_config)
    
    if not runtime_result["success"]:
        return {
            "nextflow_command": "",
            "execution_status": "failed",
            "execution_output": f"配置生成失败: {runtime_result['error']}",
            "execution_result": {"success": False, "error": runtime_result["error"]},
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
        # 资源配置摘要
        resource_summary = ""
        if resource_config:
            total_cpus = sum(config.get('cpus', 1) for config in resource_config.values())
            resource_summary = f"\n   - 资源分配: {len(resource_config)}个进程，总CPU {total_cpus}核"
        
        response_msg = f"""🎉 **RNA-seq分析执行成功！**

📋 **执行摘要:**
   - 基因组版本: {nextflow_config.get('genome_version', 'unknown')}
   - 分析工具链: {nextflow_config.get('qc_tool', 'unknown')}-{nextflow_config.get('align_tool', 'unknown')}-{nextflow_config.get('quant_tool', 'unknown')}
   - 执行时长: {execution_result.get('duration', 'unknown')}{resource_summary}
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
        "status": "analysis"  
        }