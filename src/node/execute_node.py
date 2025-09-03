import json
import subprocess
import asyncio
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional
from jinja2 import Environment, FileSystemLoader, Template
from ..state import AgentState
from ..config import get_tools_config

def _get_nextflow_template() -> Template:
    """获取 Nextflow 配置模板"""
    config = get_tools_config()
    templates_dir = config.settings.templates_dir  # 使用新的Docker兼容路径
    
    # 确保模板目录和文件存在
    templates_dir.mkdir(parents=True, exist_ok=True)
    template_file = templates_dir / "nextflow_config.j2"
    
    if not template_file.exists():
        raise FileNotFoundError(f"模板文件不存在: {template_file}")
    
    env = Environment(loader=FileSystemLoader(templates_dir))
    return env.get_template("nextflow_config.j2")

async def generate_nextflow_config(resource_config: Dict[str, Dict[str, Any]], report_dir: Optional[str] = None) -> Dict[str, Any]:
    """生成动态的nextflow.config文件，使用Jinja2模板并存放在时间戳目录中"""
    try:
        config = get_tools_config()
        
        # 使用报告目录作为配置文件存放位置
        if report_dir:
            config_dir = Path(report_dir)
        else:
            # 默认使用reports目录下的时间戳目录
            report_ts = datetime.now().strftime('%Y%m%d_%H%M%S')
            config_dir = config.reports_dir / report_ts
        
        config.path_manager.ensure_directory(config_dir)
        
        # 设置 Nextflow 报告路径（在配置目录下创建 nextflow 子目录）
        nf_reports_dir = config_dir / "nextflow"
        nf_reports_dir.mkdir(parents=True, exist_ok=True)
        report_path = nf_reports_dir / "execution_report.html"
        timeline_path = nf_reports_dir / "execution_timeline.html"
        trace_path = nf_reports_dir / "execution_trace.txt"
        
        # 准备模板变量
        template_vars = {
            "resource_config": resource_config or {},
            "report_file": report_path.as_posix(),
            "timeline_file": timeline_path.as_posix(),
            "trace_file": trace_path.as_posix(),
            # 可选的自定义参数
            "default_cpus": 1,
            "default_memory": "2 GB",
            "max_retries": 3,
            "max_errors": "-1",
            "executor_name": "local"
        }
        
        # 生成配置内容
        template = _get_nextflow_template()
        config_content = template.render(**template_vars)
        
        # 写入配置文件
        config_file = config_dir / "nextflow.config"
        with open(config_file, 'w', encoding='utf-8') as f:
            f.write(config_content)
        
        # 计算配置统计
        total_processes = len(resource_config) if resource_config else 0  
        total_cpus = sum(cfg.get('cpus', 1) for cfg in (resource_config.values() if resource_config else []))
        
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
        config = get_tools_config()
        if report_dir:
            base_dir = Path(report_dir)
        else:
            base_dir = config.reports_dir
        config.path_manager.ensure_directory(base_dir)
        
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

def build_nextflow_command(params_file_path: str, config_file_path: Optional[str] = None) -> str:
    """构建Nextflow命令（使用 params-file 传递参数，支持自定义配置文件路径）"""
    config = get_tools_config()
    
    # 优先使用传入的配置文件路径，否则使用默认路径
    if config_file_path:
        nextflow_config = config_file_path
    else:
        nextflow_config = str(config.settings.nextflow_config_path)
    
    cmd_parts = [
        "nextflow", "run", "/main.nf",
        "-c", nextflow_config,
        "-params-file", params_file_path,
        "-work-dir", "work",
        "-resume",
    ]
    return " ".join(cmd_parts)

async def execute_nextflow_pipeline(command: str) -> Dict[str, Any]:
    """执行Nextflow流水线（已优化，移除SSL重试机制）"""
    
    start_time = time.time()
    
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
        
        # 类型断言确保stdout不为None
        assert process.stdout is not None, "stdout should not be None when PIPE is specified"
        
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
        
        # 执行结果判断
        if process.returncode == 0:
            print(f"✅ Nextflow执行成功 (耗时: {duration_str})")
            return {
                "success": True,
                "output": full_output,
                "duration": duration_str,
                "return_code": process.returncode,
                "mode": "success"
            }
        else:
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
        error_msg = f"执行异常: {str(e)}"
        print(f"❌ {error_msg}")
        return {
            "success": False,
            "output": "",
            "duration": f"{duration:.1f}秒",
            "error": error_msg,
            "mode": "exception"
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
    config = get_tools_config()
    report_dir = config.reports_dir / report_ts
    config.path_manager.ensure_directory(report_dir)
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
    
    # 构建Nextflow命令（使用报告目录中的配置文件和参数文件）
    print(f"\n🔧 **构建Nextflow命令...**")
    params_file = str(report_dir / "runtime_config.json")
    config_file = config_generation_result.get("config_file")  # 获取生成的配置文件路径
    nextflow_command = build_nextflow_command(params_file, config_file_path=config_file)
    print(f"📋 命令: {nextflow_command}")
    print(f"📄 配置文件: {config_file}")
    print(f"📦 参数文件: {params_file}")
    
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
