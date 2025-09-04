from typing import Dict, Any
from ..state import AgentState

async def user_confirm_node(state: AgentState) -> Dict[str, Any]:
    """用户确认节点 - 展示配置并等待用户决策"""
    print(f"\n{'='*60}")
    print(f"🎯 **分析配置确认**")
    print(f"{'='*60}")
    
    # 展示当前配置摘要
    nextflow_config = state.nextflow_config or {}
    resource_config = state.resource_config or {}
    config_reasoning = state.config_reasoning or "系统自动生成配置"
    
    print(f"\n📋 **配置摘要:**")
    if nextflow_config:
        for key, value in nextflow_config.items():
            # 格式化显示配置项
            if key == "genome_version":
                print(f"   🧬 基因组版本: {value}")
            elif key == "species":
                print(f"   🔬 物种: {value}")
            elif key == "qc_tool":
                print(f"   🧹 质控工具: {value}")
            elif key == "align_tool":
                print(f"   🎯 比对工具: {value}")
            elif key == "quant_tool":
                print(f"   📊 定量工具: {value}")
            elif key == "sample_groups":
                print(f"   📂 样本文件: {len(value)}个样本")
                for i, sample in enumerate(value, 1):
                    sample_id = sample.get('sample_id', 'Unknown')
                    read1 = sample.get('read1', '')
                    read2 = sample.get('read2', '')
                    print(f"      {i}. {sample_id}")
                    print(f"         R1: {read1}")
                    if read2:
                        print(f"         R2: {read2}")
            elif key == "paired_end":
                end_type = "双端测序" if value else "单端测序"
                print(f"   🔄 测序类型: {end_type}")
            elif key == "run_download_genome":
                download_status = "是" if value else "否"
                print(f"   ⬇️ 下载基因组: {download_status}")
            elif key == "run_build_star_index":
                # 只有当比对工具是STAR时才显示STAR索引构建状态
                align_tool = nextflow_config.get("align_tool", "").lower()
                if align_tool == "star":
                    build_status = "是" if value else "否"
                    print(f"   🏗️ 构建STAR索引: {build_status}")
            elif key == "run_build_hisat2_index":
                # 只有当比对工具是HISAT2时才显示HISAT2索引构建状态  
                align_tool = nextflow_config.get("align_tool", "").lower()
                if align_tool == "hisat2":
                    build_status = "是" if value else "否"
                    print(f"   🏗️ 构建HISAT2索引: {build_status}")
            else:
                print(f"   ⚙️ {key}: {value}")
    else:
        print(f"   ⚠️ 无配置信息")
    
    # 展示资源配置
    if resource_config:
        print(f"\n🖥️ **资源配置:**")
        for process_name, config in resource_config.items():
            cpus = config.get('cpus', 'N/A')
            memory = config.get('memory', 'N/A')
            reasoning = config.get('reasoning', '')
            
            # 格式化进程名称显示
            display_name = {
                'prepare_star_index': '🏗️ STAR索引构建',
                'prepare_hisat2_index': '🏗️ HISAT2索引构建',
                'run_alignment': '🎯 序列比对', 
                'run_quality_control': '🧹 质控处理',
                'run_quantification': '📊 基因定量',
                'download_genome_fasta': '⬇️ FASTA下载',
                'download_genome_gtf': '⬇️ GTF下载'
            }.get(process_name, f'⚙️ {process_name}')
            
            print(f"   {display_name}: {cpus}核, {memory}")
            if reasoning:
                print(f"      💭 {reasoning}")
    else:
        print(f"\n🖥️ **资源配置:** 使用默认设置")

    # 展示 fastp 参数对比（如果存在） - 支持迭代优化显示
    try:
        fp_prev = getattr(state, 'fastp_prev_params', {}) or {}
        fp_current = getattr(state, 'fastp_current_params', {}) or {}
        fp_opt = getattr(state, 'fastp_optimized_params', {}) or {}
        fp_applied = getattr(state, 'fastp_applied_updates', {}) or {}
        fp_version = getattr(state, 'fastp_version', 1)
        fp_history = getattr(state, 'fastp_version_history', []) or []
        qc_tool = (nextflow_config.get('qc_tool') or '').lower()
        
        if qc_tool == 'fastp' and (fp_current or fp_opt):
            print(f"\n🧹 **Fastp 参数管理 [v{fp_version}]:**")
            
            # 显示版本历史摘要
            if fp_history:
                print(f"   📚 版本历史: {len(fp_history)} 个版本")
                recent_versions = fp_history[-3:] if len(fp_history) > 3 else fp_history
                for record in recent_versions:
                    v = record.get("version", "?")
                    success_rate = record.get("execution_result", {}).get("success_rate", 0)
                    param_count = len(record.get("params", {}))
                    opt_count = len(record.get("optimized_params", {}))
                    print(f"     v{v}: {param_count}参数 -> {opt_count}优化 (成功率: {success_rate:.1%})")
            
            # 仅展示与当前参数不同的优化项（使用本次实际应用的差异，避免被“已应用”导致的空显示）
            diff_opt = dict(fp_applied)

            if fp_current:
                print("   • 当前参数 (历史进化结果):")
                for k in sorted(fp_current.keys()):
                    v = fp_current[k]
                    if k in fp_applied:
                        prev_v = fp_prev.get(k, None)
                        if prev_v is None:
                            print(f"     - {k}: {v}")
                        elif prev_v != v:
                            print(f"     - {k}: {prev_v} -> {v}")
                        else:
                            print(f"     - {k}: {v}")
                    else:
                        print(f"     - {k}: {v}")
            else:
                print("   • 当前参数: 使用内置默认")

            # 仅展示与当前参数不同的优化项（分栏摘要）
            if diff_opt:
                print("   • 优化参数 (本次建议):")
                for k in sorted(diff_opt.keys()):
                    new_v = diff_opt[k]
                    old_v = fp_prev.get(k, None)
                    if old_v is not None and old_v != new_v:
                        print(f"     - {k}: {old_v} -> {new_v}")
                    else:
                        print(f"     - {k}: {new_v}")
            else:
                print("   • 优化参数: 暂无新建议")

    except Exception as _:
        # 展示失败不影响确认流程
        pass

    print(f"\n💭 **配置理由:**")
    print(f"   {config_reasoning}")
    
    print(f"\n🔄 **请选择下一步操作:**")
    print(f"   /execute_once    - ▶️ 单次执行（仅运行fastp质控）")
    print(f"   /execute_opt     - ⚡ 优化执行（运行fastp并给出组级优化建议）")
    print(f"   /modify [需求]   - 🔧 修改配置")  
    print(f"   /cancel          - ❌ 取消分析返回普通模式")
    print(f"   /quit            - 🚪 退出程序")
    print(f"{'='*60}")
    
    # 获取用户输入
    user_choice = ""  # 初始化变量避免引用错误
    try:
        user_choice = input("请输入命令: ").strip()
        
        # 处理用户输入 - 简化逻辑
        user_choice_lower = user_choice.lower()
        
        # 定义modify等价命令
        modify_prefixes = ['/modify', '/修改', '/调整']
        is_modify_command = (user_choice_lower in modify_prefixes or 
                           any(user_choice_lower.startswith(f"{prefix} ") for prefix in modify_prefixes))
        
        if user_choice_lower in ['/execute_once', '/once', '/单次执行', '/单次']:
            user_decision = "execute"
            execution_mode = 'single'
            decision_msg = "✅ 单次执行分析"
        elif user_choice_lower in ['/execute_opt', '/optimize', '/optimized', '/优化执行', '/优化']:
            user_decision = "execute"
            execution_mode = 'optimized'
            decision_msg = "⚡ 优化执行分析"
        elif user_choice_lower in ['/quit', '/exit', '/退出', '/bye']:
            user_decision = "quit"
            decision_msg = "🚪 退出程序"
        elif is_modify_command:
            user_decision = "modify"
            decision_msg = "🔧 修改配置"
            
            # 处理modify等价命令中的新需求 - 优雅的参数提取
            modify_content = ""
            for prefix in modify_prefixes:
                if user_choice_lower.startswith(prefix):
                    modify_content = user_choice_lower.replace(prefix, '', 1).strip()
                    break
            
            if modify_content:
                new_user_requirements = {"raw_input": modify_content}
            else:
                new_user_requirements = {}
        elif user_choice_lower in ['/cancel', '/取消']:
            user_decision = "cancel"
            decision_msg = "❌ 取消分析"
        else:
            # 无效输入，提示用户重新选择
            print(f"❌ 无效输入: {user_choice}")
            print(f"请选择有效的命令: /execute_once, /execute_opt, /modify, /cancel, /quit")
            # 递归调用自己，重新获取用户输入
            return await user_confirm_node(state)
        
        print(f"🎯 {decision_msg}")
        
    except KeyboardInterrupt:
        print(f"\n⚠️ 用户中断，取消分析")
        user_choice = "/cancel"  # 设置默认值避免引用错误
        user_decision = "cancel"
        decision_msg = "❌ 用户中断取消"
    except Exception as e:
        print(f"❌ 输入处理错误: {e}")
        user_choice = "/cancel"  # 设置默认值避免引用错误
        user_decision = "cancel"
        decision_msg = "❌ 输入错误取消"
    
    confirmation_message = f"""🎯 分析配置已确认！

📋 配置项: {len(nextflow_config)} 个参数已设置
💭 决策理由: {config_reasoning}
🎯 用户选择: {decision_msg}

准备进入下一阶段..."""
    
    return {
        # 从prepare_node继承并传递
        "nextflow_config": nextflow_config,
        "resource_config": resource_config,
        "config_reasoning": config_reasoning,
        
        # 当前节点输出
        "confirmation_message": confirmation_message,
        "user_decision": user_decision,
        "response": decision_msg,
        "status": user_decision,
        "execution_mode": locals().get('execution_mode', getattr(state, 'execution_mode', 'single')),
        
        # 重新修改时设置modify需求，保持初始user_requirements不变
        "user_requirements": getattr(state, 'user_requirements', {}),  # 保持初始需求
        "modify_requirements": new_user_requirements if 'new_user_requirements' in locals() else {},  # modify需求
        
        # 保存用户选择用于后续处理
        "messages": [{"role": "user", "content": user_choice}]
    }
