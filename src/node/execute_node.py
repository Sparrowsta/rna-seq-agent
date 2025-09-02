import json
import subprocess
import asyncio
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional
from ..state import AgentState

async def generate_nextflow_config(resource_config: Dict[str, Dict[str, Any]], report_dir: Optional[str] = None) -> Dict[str, Any]:
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

        # 添加LLM智能分配的资源配置（仅使用 prepare 提供的进程名；无默认 withName 片段）
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

        # 添加执行器和报告配置（报告目录指向 reports/<ts>/nextflow/ 或默认 results/nextflow/）
        nf_reports_dir = Path(report_dir) / "nextflow" if report_dir else Path("results/nextflow")
        report_path = nf_reports_dir / "execution_report.html"
        timeline_path = nf_reports_dir / "execution_timeline.html"
        trace_path = nf_reports_dir / "execution_trace.txt"

        config_content += f"""}}

// 执行配置
executor {{
    name = 'local'
}}

// 报告配置
report {{
    enabled = true
    file = '{report_path.as_posix()}'
    overwrite = true
}}

timeline {{
    enabled = true
    file = '{timeline_path.as_posix()}'
    overwrite = true
}}

trace {{
    enabled = true
    file = '{trace_path.as_posix()}'
    overwrite = true
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,read_bytes,write_bytes'
}}
"""

        # 写入配置文件
        config_file = config_dir / "nextflow.config"
        with open(config_file, 'w', encoding='utf-8') as f:
            f.write(config_content)
        
        # 计算配置统计（无默认 withName，空则为 0）
        total_processes = len(resource_config) if resource_config else 0
        total_cpus = sum(cfg.get('cpus', 1) for cfg in resource_config.values()) if resource_config else 0
        
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

async def generate_runtime_config(nextflow_config: Dict[str, Any], resource_config: Optional[Dict[str, Dict[str, Any]]] = None, report_dir: Optional[str] = None) -> Dict[str, Any]:
    """生成运行时配置文件"""
    try:
        base_dir = Path(report_dir) if report_dir else Path("reports")
        base_dir.mkdir(parents=True, exist_ok=True)
        
        # 明确处理resource_config的None情况
        if resource_config is None:
            resource_config = {}
        
        # 创建运行时配置（以 Nextflow params-file 直读的扁平键为准）
        flat_params: Dict[str, Any] = dict(nextflow_config or {})
        runtime_config = {
            **flat_params,
            "__meta": {
                "timestamp": datetime.now().isoformat(),
                "analysis_id": f"rna_seq_{int(time.time())}",
            },
            # 保留资源配置供归档与追溯（Nextflow 会忽略该未知参数）
            "resource_config": resource_config,
        }
        
        # 保存配置文件
        config_file = base_dir / "runtime_config.json"
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(runtime_config, f, indent=2, ensure_ascii=False)
        
        print(f"✅ 配置文件已保存: {config_file}")
        print(f"📦 传递给Nextflow的参数: {list(flat_params.keys())}")
        return {"success": True, "config_file": str(config_file)}
        
    except Exception as e:
        print(f"❌ 配置生成失败: {e}")
        return {"success": False, "error": str(e)}

def build_nextflow_command(nextflow_config: Dict[str, Any], params_file_path: str) -> str:
    """构建Nextflow命令（使用 params-file 传递参数）"""
    cmd_parts = [
        "nextflow", "run", "/main.nf",
        "-c", "/config/nextflow.config",
        "-params-file", params_file_path,
        "-work-dir", "work",
        "-resume",
    ]
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
    
    # 生成报告时间戳与目录
    report_ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_dir = Path("reports") / report_ts
    report_dir.mkdir(parents=True, exist_ok=True)
    # 确保 Nextflow 报告子目录存在，避免部分环境下不自动创建
    nf_reports_dir = report_dir / "nextflow"
    nf_reports_dir.mkdir(parents=True, exist_ok=True)
    print(f"📁 报告目录: {report_dir}")

    # 生成动态的nextflow.config文件（报告目录传入以定向 Nextflow 报告输出）
    print(f"\n⚙️ **生成Nextflow配置文件...**")
    config_generation_result = await generate_nextflow_config(resource_config, report_dir=str(report_dir))
    
    if not config_generation_result["success"]:
        return {
            "nextflow_command": "",
            "execution_status": "failed",
            "execution_output": f"Nextflow配置生成失败: {config_generation_result['error']}",
            "execution_result": {"success": False, "error": config_generation_result["error"]},
            "response": "分析执行失败：Nextflow配置文件生成错误",
            "status": "failed"
        }
    
    # 生成运行时配置文件（写入报告目录）
    print(f"\n📝 **生成运行时配置...**")
    runtime_result = await generate_runtime_config(nextflow_config, resource_config, report_dir=str(report_dir))
    
    if not runtime_result["success"]:
        return {
            "nextflow_command": "",
            "execution_status": "failed",
            "execution_output": f"配置生成失败: {runtime_result['error']}",
            "execution_result": {"success": False, "error": runtime_result["error"]},
            "response": "分析执行失败：配置生成错误",
            "status": "failed"
        }
    
    # 构建Nextflow命令（使用报告目录中的 params-file）
    print(f"\n🔧 **构建Nextflow命令...**")
    params_file = str(report_dir / "runtime_config.json")
    nextflow_command = build_nextflow_command(nextflow_config, params_file_path=params_file)
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
        "report_dir": str(report_dir),
        "report_ts": report_ts,
        "response": response_msg,
        "status": "analysis"  
        }
