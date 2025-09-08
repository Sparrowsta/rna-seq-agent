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

    # =============== 三层对比展示 ===============
    def _flatten(d: Dict[str, Any], parent: str = "", sep: str = ".") -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        if isinstance(d, dict):
            for k, v in d.items():
                key = f"{parent}{sep}{k}" if parent else str(k)
                if isinstance(v, dict):
                    out.update(_flatten(v, key, sep))
                else:
                    out[key] = v
        return out

    # 1) Base 快照（首次进入确认页时建立基线）
    base_nextflow_config = getattr(state, 'prepare_defaults_nextflow_config', None) or {}
    base_resource_config = getattr(state, 'prepare_defaults_resource_config', None) or {}
    base_fastp_params = getattr(state, 'prepare_defaults_fastp_params', None) or {}

    if not base_nextflow_config:
        base_nextflow_config = dict(nextflow_config)
    if not base_resource_config:
        base_resource_config = dict(resource_config)
    if not base_fastp_params:
        base_fastp_params = dict(getattr(state, 'fastp_params', {}) or {})

    # 2) Effective 当前生效值
    effective_nextflow_config = dict(nextflow_config)
    effective_resource_config = dict(resource_config)
    effective_fastp_params = dict(getattr(state, 'fastp_params', {}) or {})

    # 3) Mods（根据 Base 与 Effective 的差异推断，Nextflow/Resource 仅来源于用户修改）
    ignore_nf_keys = {"results_dir", "results_timestamp", "base_data_path", "validated_work_dir", "sample_groups"}
    flattened_base_nextflow = {k: v for k, v in _flatten(base_nextflow_config).items() if k.split('.')[0] not in ignore_nf_keys and k not in ignore_nf_keys}
    flattened_effective_nextflow = {k: v for k, v in _flatten(effective_nextflow_config).items() if k.split('.')[0] not in ignore_nf_keys and k not in ignore_nf_keys}
    nextflow_modifications = {k: (flattened_base_nextflow.get(k), flattened_effective_nextflow.get(k)) for k in flattened_effective_nextflow.keys() | flattened_base_nextflow.keys() if flattened_base_nextflow.get(k) != flattened_effective_nextflow.get(k)}

    # 资源配置差异（逐进程关注 cpus/memory）
    resource_modifications: Dict[str, Dict[str, Any]] = {}
    for process_name_key in set(base_resource_config.keys()) | set(effective_resource_config.keys()):
        base_process_config = base_resource_config.get(process_name_key, {}) or {}
        effective_process_config = effective_resource_config.get(process_name_key, {}) or {}
        process_diff: Dict[str, Any] = {}
        for key in {"cpus", "memory"}:
            if base_process_config.get(key) != effective_process_config.get(key):
                process_diff[key] = (base_process_config.get(key), effective_process_config.get(key))
        if process_diff:
            resource_modifications[process_name_key] = process_diff

    # 4) Opt（FastP 优化建议）
    fastp_history = getattr(state, 'fastp_params_history', []) or []
    last_applied = fastp_history[-1].get('optimization_applied') if fastp_history else {}
    fastp_opt = getattr(state, 'fastp_optimized_suggestions', {}) or {}

    # =============== 展示 ===============
    # Nextflow 配置 Mods 差异
    if nextflow_modifications:
        print(f"\n🧭 **配置对比（Nextflow）**")
        for config_key in sorted(nextflow_modifications.keys()):
            old_value, new_value = nextflow_modifications[config_key]
            print(f"   - {config_key}: {old_value} -> {new_value}")
    # 资源配置 Mods 差异
    if resource_modifications:
        print(f"\n🧮 **配置对比（资源）**")
        for process_name_key, process_diff in resource_modifications.items():
            for resource_key, (old_value, new_value) in process_diff.items():
                print(f"   - {process_name_key}.{resource_key}: {old_value} -> {new_value}")

    # FastP 三层对比
    try:
        qc_tool = (nextflow_config.get('qc_tool') or '').lower()
        if qc_tool == 'fastp' and effective_fastp_params:
            print(f"\n🧹 **FastP 参数（三层对比）**")
            # Effective（当前）
            print(f"   📋 当前（Effective）:")
            for param_key in sorted(effective_fastp_params.keys()):
                print(f"     - {param_key}: {effective_fastp_params[param_key]}")

            # Mods：与 Base 对比的差异，排除已由优化应用的键
            flattened_base_fastp = _flatten(base_fastp_params)
            flattened_effective_fastp = _flatten(effective_fastp_params)
            applied_keys = set((last_applied or {}).keys())

            # 优先显示最近一次“用户修改”的键；若没有，则按差异回退
            user_modified_keys: list[str] = []
            try:
                modification_history = getattr(state, 'modification_history', []) or []
                if modification_history:
                    last_record = modification_history[-1] or {}
                    last_fastp_changes = (last_record.get('changes') or {}).get('fastp') or {}
                    if isinstance(last_fastp_changes, dict):
                        user_modified_keys = list(last_fastp_changes.keys())
            except Exception:
                user_modified_keys = []

            if user_modified_keys:
                print(f"\n   ✏️ 用户修改（Mods）:")
                for param_key in sorted(user_modified_keys):
                    print(f"     - {param_key}: {flattened_base_fastp.get(param_key)} -> {flattened_effective_fastp.get(param_key)}")
            else:
                modified_keys = [
                    k for k in flattened_effective_fastp.keys()
                    if flattened_base_fastp.get(k) != flattened_effective_fastp.get(k) and k not in applied_keys
                ]
                if modified_keys:
                    print(f"\n   ✏️ 用户修改（Mods）:")
                    for param_key in sorted(modified_keys):
                        print(f"     - {param_key}: {flattened_base_fastp.get(param_key)} -> {flattened_effective_fastp.get(param_key)}")
                else:
                    print(f"\n   ✏️ 用户修改（Mods）: 无")

            # Opt：优化建议（区分已应用/待应用）
            if fastp_opt:
                print(f"\n   ⚙️ 优化建议（Opt）:")
                for param_key in sorted(fastp_opt.keys()):
                    suggestion_value = fastp_opt[param_key]
                    flag = "[applied]" if param_key in applied_keys and flattened_effective_fastp.get(param_key) == suggestion_value else "[pending]"
                    current_value = flattened_effective_fastp.get(param_key)
                    if flag == "[applied]":
                        print(f"     - {param_key}: {current_value} (已应用 {suggestion_value}) {flag}")
                    else:
                        print(f"     - {param_key}: {current_value} -> {suggestion_value} {flag}")
            else:
                print(f"\n   ⚙️ 优化建议（Opt）: 无")
                
            # 显示最新的优化建议（如果有的话）
            latest_optimization_suggestions = getattr(state, 'fastp_optimization_suggestions', {})
            latest_optimization_reasoning = getattr(state, 'fastp_optimization_reasoning', '')
            
            if latest_optimization_suggestions:
                print(f"\n   🔍 **最新优化分析**:")
                print(f"     💡 {latest_optimization_reasoning}")
                print(f"     📊 具体建议:")
                for param_key, suggested_value in latest_optimization_suggestions.items():
                    current_val = flattened_effective_fastp.get(param_key, "未设置")
                    print(f"       - {param_key}: {current_val} → {suggested_value}")
    except Exception as _:
        pass

    print(f"\n💭 **配置理由:**")
    print(f"   {config_reasoning}")
    
    print(f"\n🔄 **请选择下一步操作:**")
    
    # 根据执行进度显示不同的选项
    completed_steps = getattr(state, 'completed_steps', [])
    current_step = getattr(state, 'current_step', '')
    
    if completed_steps:
        print(f"   📊 **执行进度**: {' -> '.join(completed_steps)}")
        if current_step:
            print(f"   🔄 **当前步骤**: {current_step}")
        
        # 根据进度显示继续选项
        if "featurecounts" in completed_steps:
            print(f"   /continue        - ➡️ 继续到综合分析")
        elif "star" in completed_steps:
            print(f"   /continue        - ➡️ 继续到FeatureCounts定量")
        elif "fastp" in completed_steps:
            print(f"   /continue        - ➡️ 继续到STAR比对")
        
        print(f"   /restart         - 🔄 重新开始完整流水线")
    
    # 检查是否有优化建议可以应用
    latest_optimization_suggestions = getattr(state, 'fastp_optimization_suggestions', {})
    if latest_optimization_suggestions:
        print(f"   /apply_opt       - ✨ 应用最新优化建议")
    
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
        elif user_choice_lower in ['/continue', '/继续']:
            # 根据当前进度决定下一步 - 只有有进度时才允许continue
            if not completed_steps:
                # 没有任何进度，不能continue
                print(f"❌ 无执行进度，请先选择 /execute_once 或 /execute_opt 开始分析")
                return await user_confirm_node(state)
            elif "featurecounts" in completed_steps:
                user_decision = "continue_analysis"
                decision_msg = "➡️ 继续到综合分析"
            elif "star" in completed_steps:
                user_decision = "continue_featurecounts"
                decision_msg = "➡️ 继续到FeatureCounts定量"
            elif "fastp" in completed_steps:
                user_decision = "continue_star"
                decision_msg = "➡️ 继续到STAR比对"
        elif user_choice_lower in ['/restart', '/重启', '/重新开始']:
            user_decision = "execute"
            execution_mode = 'single'
            decision_msg = "🔄 重新开始完整流水线"
            # 清空进度信息
            completed_steps = []
            current_step = ""
        elif user_choice_lower in ['/apply_opt', '/应用优化', '/优化应用']:
            # 应用最新的优化建议
            latest_optimization_suggestions = getattr(state, 'fastp_optimization_suggestions', {})
            if latest_optimization_suggestions:
                user_decision = "apply_optimization"
                decision_msg = f"✨ 应用优化建议（{len(latest_optimization_suggestions)}项参数）"
            else:
                print(f"❌ 当前没有可用的优化建议")
                return await user_confirm_node(state)
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
        
        # 进度信息
        "completed_steps": locals().get('completed_steps', getattr(state, 'completed_steps', [])),
        "current_step": locals().get('current_step', getattr(state, 'current_step', '')),
        
        # 重新修改时设置modify需求，保持初始user_requirements不变
        "user_requirements": getattr(state, 'user_requirements', {}),  # 保持初始需求
        "modify_requirements": new_user_requirements if 'new_user_requirements' in locals() else {},  # modify需求
        
        # 保存用户选择用于后续处理
        "messages": [{"role": "user", "content": user_choice}],
        # 持久化 Base 快照（便于后续对比）
        "prepare_defaults_nextflow_config": base_nextflow_config,
        "prepare_defaults_resource_config": base_resource_config,
        "prepare_defaults_fastp_params": base_fastp_params
    }
