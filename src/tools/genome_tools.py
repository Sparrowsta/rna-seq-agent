"""
RNA-seq智能分析助手 - 基因组管理工具

包含：
- add_genome_config: 添加基因组配置
- download_genome_assets: 下载基因组文件
- build_star_index: 构建STAR索引
- build_hisat2_index: 构建HISAT2索引
"""

import gzip
import json
import shutil
import subprocess
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

# 使用官方工具装饰器
from langchain_core.tools import tool

# 导入配置模块
from ..config import get_tools_config
from ..logging_bootstrap import get_logger

logger = get_logger("rna.tools.genome")
@tool
def add_genome_config(genome_info: Dict[str, Any]) -> Dict[str, Any]:
    """添加基因组配置到 genomes.json（不在工具内调用LLM）。

    期望输入（由 Normal 节点的 LLM 已解析好并传入）：
    {
      "genome_id": "如 hg38/mm39/...",
      "species": "如 human/mouse/...",
      "version": "通常与 genome_id 相同",
      "fasta_url": "FASTA 完整下载 URL",
      "gtf_url": "GTF 完整下载 URL",
      // 可选：当需要自定义时可传入以下两项；否则自动按规范路径生成
      "fasta_path": "genomes/<species>/<version>/<version>.fa",
      "gtf_path": "genomes/<species>/<version>/<version>.gtf"
    }
    """
    try:
        # 兼容：若传入字符串，尝试按 JSON 解析
        if isinstance(genome_info, str):
            try:
                genome_info = json.loads(genome_info)
            except Exception:
                return {"success": False, "error": "genome_info 需为对象或可解析为JSON的字符串"}

        if not genome_info or not isinstance(genome_info, dict):
            return {"success": False, "error": "参数 genome_info 不能为空，且需为对象"}

        # 验证必需字段
        required_fields = ["genome_id", "species", "version", "fasta_url", "gtf_url"]
        missing_fields = [f for f in required_fields if not genome_info.get(f)]
        if missing_fields:
            return {"success": False, "error": f"缺少字段：{missing_fields}"}

        genome_id = str(genome_info["genome_id"]).strip()
        species = str(genome_info["species"]).strip()
        version = str(genome_info["version"]).strip()
        if not genome_id or not species or not version:
            return {"success": False, "error": "genome_id/species/version 不能为空"}

        # 规范化本地存放路径（允许外部覆盖）
        fasta_path_raw = genome_info.get("fasta_path") or f"genomes/{species}/{version}/{version}.fa"
        gtf_path_raw = genome_info.get("gtf_path") or f"genomes/{species}/{version}/{version}.gtf"
        
        fasta_path = fasta_path_raw
        gtf_path = gtf_path_raw

        new_genome_config = {
            "species": species,
            "version": version,
            "fasta_path": fasta_path,
            "gtf_path": gtf_path,
            "fasta_url": genome_info["fasta_url"],
            "gtf_url": genome_info["gtf_url"]
        }

        # 读取现有配置并添加
        config = get_tools_config()
        genomes_file = config.genomes_config_path

        # 如文件不存在，初始化为空对象（尽量自愈）
        if genomes_file.exists():
            with open(genomes_file, "r", encoding="utf-8") as f:
                try:
                    genomes_data = json.load(f) or {}
                except json.JSONDecodeError:
                    genomes_data = {}
        else:
            genomes_file.parent.mkdir(parents=True, exist_ok=True)
            genomes_data = {}

        # 检查重复
        if genome_id in genomes_data:
            logger.warning(f"重复的基因组ID：{genome_id}")
            return {
                "success": False,
                "error": f"基因组 {genome_id} 已存在",
                "existing_config": genomes_data[genome_id]
            }

        # 保存配置
        genomes_data[genome_id] = new_genome_config

        with open(genomes_file, "w", encoding="utf-8") as f:
            json.dump(genomes_data, f, indent=2, ensure_ascii=False)

        logger.info(f"已添加基因组配置：{genome_id} -> {genomes_file}")
        return {
            "success": True,
            "genome_id": genome_id,
            "config": new_genome_config,
            "config_path": str(genomes_file),
            "message": f"成功添加基因组配置：{genome_id}"
        }

    except Exception as e:
        logger.error(f"添加基因组配置失败：{e}")
        return {"success": False, "error": f"添加基因组配置时出错：{str(e)}"}


def _download_genome_file(
    url: str,
    destination: Path,
    *,
    timeout: int,
    label: str,
    timeout_message: str,
) -> Tuple[bool, Optional[str], Optional[Path]]:
    """从远程下载基因组文件到 tmp 目录，成功后移动到目标位置。

    返回：
        (success, error_message, tmp_file_path)
        - success: 是否成功
        - error_message: 若失败则包含错误详情
        - tmp_file_path: 实际下载文件路径，便于排查
    """
    tools_config = get_tools_config()
    tmp_dir = Path(tools_config.settings.data_dir) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    token = uuid.uuid4().hex
    download_tmp = tmp_dir / f"{destination.name}.{token}"
    prepared_tmp = download_tmp
    needs_decompress = url.endswith(".gz") and not destination.name.endswith(".gz")

    if needs_decompress:
        download_tmp = tmp_dir / f"{destination.name}.{token}.gz"
        prepared_tmp = tmp_dir / f"{destination.name}.{token}.ready"

    try:
        result = subprocess.run(
            [
                "curl",
                "-L",
                "-sS",
                "--fail",
                "--retry",
                "2",
                "--connect-timeout",
                "30",
                "--max-time",
                str(timeout),
                "-o",
                str(download_tmp),
                url,
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        if result.returncode != 0:
            error_detail = result.stderr.strip() or result.stdout.strip() or f"curl exit code {result.returncode}"
            return False, f"{label}下载失败: {error_detail} (临时文件: {download_tmp})", download_tmp

        if needs_decompress:
            try:
                with gzip.open(download_tmp, "rb") as src, open(prepared_tmp, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                logger.info(f"{label}下载完成并解压: {prepared_tmp}")
            except Exception as exc:  # noqa: BLE001
                return False, f"{label}解压失败: {exc} (临时文件: {download_tmp})", download_tmp
        else:
            prepared_tmp = download_tmp
            logger.info(f"{label}下载完成: {prepared_tmp}")

        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(prepared_tmp), str(destination))
            logger.info(f"{label}文件已移动到目标路径: {destination}")
        except Exception as exc:  # noqa: BLE001
            return False, f"{label}移动到目标路径失败: {exc} (临时文件: {prepared_tmp})", prepared_tmp

        # 若需要解压，保留原始压缩文件供复查；否则已被移动
        return True, None, download_tmp if needs_decompress else None

    except subprocess.TimeoutExpired:
        return False, f"{timeout_message} (临时文件: {download_tmp})", download_tmp
    except Exception as exc:  # noqa: BLE001
        return False, f"{label}下载异常: {exc} (临时文件: {download_tmp})", download_tmp


@tool
def download_genome_assets(
    genome_id: str,
    download_fasta: bool = True,
    download_gtf: bool = True
) -> Dict[str, Any]:
    """下载指定基因组的FASTA和/或GTF文件

    Args:
        genome_id: 基因组标识符，如'hg38', 'mm39'等
        download_fasta: 是否下载FASTA文件，默认True
        download_gtf: 是否下载GTF文件，默认True

    Returns:
        下载结果字典，包含状态、文件路径等信息
    """
    try:
        # 参数验证：至少要选择一种文件类型
        if not download_fasta and not download_gtf:
            return {
                "success": False,
                "error": "至少需要选择下载一种文件类型（FASTA或GTF）"
            }

        config = get_tools_config()

        # 获取基因组配置
        genomes_config_path = config.genomes_config_path
        if not genomes_config_path.exists():
            return {
                "success": False,
                "error": f"基因组配置文件不存在: {genomes_config_path}"
            }
        
        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        
        if genome_id not in genomes_config:
            return {
                "success": False,
                "error": f"未找到基因组配置: {genome_id}",
                "available_genomes": list(genomes_config.keys())
            }
        
        genome_info = genomes_config[genome_id]
        fasta_url = genome_info.get("fasta_url")
        gtf_url = genome_info.get("gtf_url")
        fasta_path = genome_info.get("fasta_path")
        gtf_path = genome_info.get("gtf_path")
        
        if not all([fasta_url, gtf_url, fasta_path, gtf_path]):
            return {
                "success": False,
                "error": f"基因组配置不完整: {genome_id}",
                "config": genome_info
            }
        
        # 转换为绝对路径
        fasta_full_path = config.settings.data_dir / fasta_path
        gtf_full_path = config.settings.data_dir / gtf_path
        
        # 创建目录
        fasta_full_path.parent.mkdir(parents=True, exist_ok=True)
        gtf_full_path.parent.mkdir(parents=True, exist_ok=True)
        
        start_time = time.time()
        downloaded_files = []

        # 下载FASTA文件
        if download_fasta:
            if not fasta_full_path.exists():
                logger.info(f"下载FASTA: {fasta_url} -> {fasta_full_path}")
                success, error_msg, tmp_file = _download_genome_file(
                    fasta_url,
                    fasta_full_path,
                    timeout=3600,
                    label="FASTA",
                    timeout_message="FASTA下载超时（1小时）"
                )
                if not success:
                    return {
                        "success": False,
                        "error": error_msg,
                        "tmp_file": str(tmp_file) if tmp_file else None,
                        "execution_time": time.time() - start_time
                    }
                downloaded_files.append("fasta")
            else:
                logger.info(f"FASTA文件已存在: {fasta_full_path}")
        else:
            logger.info(f"跳过FASTA文件下载（用户选择）")

        # 下载GTF文件
        if download_gtf:
            if not gtf_full_path.exists():
                logger.info(f"下载GTF: {gtf_url} -> {gtf_full_path}")
                success, error_msg, tmp_file = _download_genome_file(
                    gtf_url,
                    gtf_full_path,
                    timeout=1800,
                    label="GTF",
                    timeout_message="GTF下载超时（30分钟）"
                )
                if not success:
                    return {
                        "success": False,
                        "error": error_msg,
                        "tmp_file": str(tmp_file) if tmp_file else None,
                        "execution_time": time.time() - start_time
                    }
                downloaded_files.append("gtf")
            else:
                logger.info(f"GTF文件已存在: {gtf_full_path}")
        else:
            logger.info(f"跳过GTF文件下载（用户选择）")

        execution_time = time.time() - start_time

        # 构建文件信息字典（只包含请求的文件）
        files_dict = {}
        skipped_by_user = []

        if download_fasta:
            files_dict["fasta"] = str(fasta_full_path)
        else:
            skipped_by_user.append("fasta")

        if download_gtf:
            files_dict["gtf"] = str(gtf_full_path)
        else:
            skipped_by_user.append("gtf")

        # 构建返回消息
        requested_types = []
        if download_fasta:
            requested_types.append("FASTA")
        if download_gtf:
            requested_types.append("GTF")
        message = f"基因组{'/'.join(requested_types)}处理完成: {genome_id}"

        return {
            "success": True,
            "genome_id": genome_id,
            "requested_files": {
                "fasta": download_fasta,
                "gtf": download_gtf
            },
            "downloaded_files": downloaded_files,
            "skipped_by_user": skipped_by_user,
            "files": files_dict,
            "tmp_file": None,
            "execution_time": execution_time,
            "message": message
        }
    
    except Exception as e:
        logger.error(f"基因组下载异常: {e}")
        return {
            "success": False,
            "error": f"基因组下载失败: {str(e)}",
            "execution_time": 0
        }


@tool
def build_star_index(
    genome_id: str,
    sjdb_overhang: Optional[int] = None,
    runThreadN: Optional[int] = None,
    force_rebuild: bool = False,
    results_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """构建STAR基因组索引（使用 build_index.nf）
    
    主要在RNA-seq分析流程中自动调用，当系统检测到缺少STAR索引时触发。
    参数文件统一保存到results/star目录，与其他STAR相关文件放在一起。

    Args:
        genome_id: 基因组标识，用于定位FASTA和GTF文件
        sjdb_overhang: 剪接位点overhang（默认用 DEFAULT_STAR_PARAMS 或 100）
        runThreadN: 线程数（默认用 DEFAULT_STAR_PARAMS 或 4）
        force_rebuild: 若索引目录已存在是否强制重建
        results_dir: 分析结果目录，参数文件保存到results/star/子目录

    Returns:
        Dict: 执行结果（包含 index_dir/stdout/stderr/skipped 等）
    """
    try:
        from ..config.default_tool_params import DEFAULT_STAR_PARAMS
        tools_config = get_tools_config()

        # 读取基因组配置
        genomes_config_path = tools_config.genomes_config_path
        if not genomes_config_path.exists():
            return {"success": False, "error": f"基因组配置文件不存在: {genomes_config_path}"}

        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genomes_config = json.load(f)
        if genome_id not in genomes_config:
            return {"success": False, "error": f"基因组 '{genome_id}' 在配置文件中不存在"}

        cfg = genomes_config[genome_id]
        fasta_rel = cfg.get("fasta_path")
        gtf_rel = cfg.get("gtf_path")
        if not fasta_rel or not gtf_rel:
            return {"success": False, "error": f"基因组 '{genome_id}' 缺少 fasta_path 或 gtf_path"}

        data_dir = tools_config.settings.data_dir
        fasta_path = data_dir / fasta_rel
        gtf_path = data_dir / gtf_rel
        
        # 检查文件是否存在，如果不存在则尝试下载
        files_missing = []
        files_invalid = []
        
        if not fasta_path.exists():
            files_missing.append(f"FASTA: {fasta_path}")
        else:
            # 检查FASTA文件有效性
            try:
                file_size = fasta_path.stat().st_size
                if file_size == 0:
                    files_invalid.append(f"FASTA文件为空: {fasta_path}")
                else:
                    # 检查文件格式（前几个字符）
                    with open(fasta_path, 'rb') as f:
                        first_bytes = f.read(10)
                    if not first_bytes.startswith(b'>'):
                        files_invalid.append(f"FASTA文件格式无效: {fasta_path} (首字符: {first_bytes[:5]})")
            except Exception as e:
                files_invalid.append(f"FASTA文件检查失败: {fasta_path} ({e})")
        
        if not gtf_path.exists():
            files_missing.append(f"GTF: {gtf_path}")
        else:
            # 检查GTF文件有效性
            try:
                file_size = gtf_path.stat().st_size
                if file_size == 0:
                    files_invalid.append(f"GTF文件为空: {gtf_path}")
            except Exception as e:
                files_invalid.append(f"GTF文件检查失败: {gtf_path} ({e})")
            
        if files_missing or files_invalid:
            logger.info(f"检测到问题文件，尝试重新下载: missing={files_missing}, invalid={files_invalid}")
            
            # 删除无效文件
            for invalid_file in files_invalid:
                try:
                    if "FASTA" in invalid_file and fasta_path.exists():
                        fasta_path.unlink()
                        logger.info(f"删除无效FASTA文件: {fasta_path}")
                    elif "GTF" in invalid_file and gtf_path.exists():
                        gtf_path.unlink()
                        logger.info(f"删除无效GTF文件: {gtf_path}")
                except Exception as e:
                    logger.warning(f"删除无效文件失败: {e}")
            
            download_result = download_genome_assets(genome_id)
            
            if not download_result.get("success", False):
                return {
                    "success": False,
                    "error": f"重新下载基因组文件失败: {download_result.get('error', '未知错误')}"
                }
            
            # 等待并再次检查文件
            time.sleep(2)  # 等待文件系统同步
            
            # 重新检查文件是否存在且有效
            if not fasta_path.exists():
                return {"success": False, "error": f"重新下载后FASTA文件仍不存在: {fasta_path}"}
            if not gtf_path.exists():
                return {"success": False, "error": f"重新下载后GTF文件仍不存在: {gtf_path}"}
            
            # 再次验证文件有效性
            try:
                fasta_size = fasta_path.stat().st_size
                if fasta_size == 0:
                    return {"success": False, "error": f"重新下载的FASTA文件仍为空: {fasta_path}"}
                
                with open(fasta_path, 'rb') as f:
                    first_bytes = f.read(10)
                if not first_bytes.startswith(b'>'):
                    return {"success": False, "error": f"重新下载的FASTA文件格式仍无效: {fasta_path}"}
                
                gtf_size = gtf_path.stat().st_size
                if gtf_size == 0:
                    return {"success": False, "error": f"重新下载的GTF文件仍为空: {gtf_path}"}
                    
            except Exception as e:
                return {"success": False, "error": f"验证重新下载的文件失败: {e}"}
            
            logger.info(f"文件重新下载并验证成功: FASTA={fasta_size}字节, GTF={gtf_size}字节")

        # 目标索引目录
        index_dir = tools_config.get_star_index_dir(fasta_path)

        # 存在即跳过
        if index_dir.exists() and not force_rebuild:
            key_files = ["SA", "SAindex", "Genome"]
            if all((index_dir / f).exists() for f in key_files):
                logger.info(f"STAR索引已存在，跳过构建: {index_dir}")
                return {
                    "success": True,
                    "index_dir": str(index_dir),
                    "skipped": True,
                    "message": "STAR索引已存在，跳过构建",
                }

        # 构造 Nextflow 参数
        sjdb = (
            int(sjdb_overhang)
            if sjdb_overhang is not None
            else int(DEFAULT_STAR_PARAMS.get("sjdbOverhang") or 100)
        )
        threads = (
            int(runThreadN)
            if runThreadN is not None
            else int(DEFAULT_STAR_PARAMS.get("runThreadN", 4))
        )

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        temp_dir = Path(tools_config.settings.data_dir) / "tmp" / f"nextflow_star_index_{timestamp}"
        temp_dir.mkdir(parents=True, exist_ok=True)

        nf_params = {
            "genome_fasta": str(fasta_path),
            "genome_gtf": str(gtf_path),
            "star_index_dir": str(index_dir),
            "sjdb_overhang": sjdb,
            "runThreadN": threads,
            "limitGenomeGenerateRAM": 32000000000,
        }

        # 创建参数文件 - 统一保存到results/star目录（与star_params在一起）
        if results_dir:
            results_path = Path(results_dir)
            star_subdir = results_path / "star"
            star_subdir.mkdir(parents=True, exist_ok=True)
            params_file = star_subdir / "build_index_params.json"
        else:
            # 罕见的独立使用情况，仍然保存到基因组star目录
            star_dir = index_dir.parent / "star"
            star_dir.mkdir(parents=True, exist_ok=True)
            params_file = star_dir / f"build_index_params_{timestamp}.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 定位 build_index.nf
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "build_index.nf"
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 build_index.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "build_index.nf")],
            }

        # 运行 Nextflow
        # 统一 Nextflow 工作目录到 /data/work
        work_root = tools_config.settings.data_dir / "work"
        work_dir = work_root / f"star_index_{timestamp}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            "nextflow",
            "run",
            str(nextflow_script),
            "-params-file",
            str(params_file),
            "-work-dir",
            str(work_dir),
        ]

        logger.info("构建STAR索引 (Nextflow)")
        logger.info(f"执行命令: {' '.join(cmd)}")
        logger.info(f"参数文件: {params_file}")
        logger.info(f"目标目录: {index_dir}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,
            cwd=tools_config.settings.project_root,
        )

        payload = {
            "success": result.returncode == 0,
            "index_dir": str(index_dir),
            "skipped": False,
            "params_file": str(params_file),
        }
        try:
            if get_tools_config().settings.debug_mode:
                payload.update({
                    "stderr": result.stderr,
                    "cmd": " ".join(cmd)
                })
        except Exception:
            pass
        # 记录结果日志
        try:
            if payload["success"]:
                logger.info(f"STAR索引构建完成：{index_dir}")
            else:
                logger.warning(f"STAR索引构建失败：rc={result.returncode} stderr={(result.stderr or '')[:400]}")
        except Exception:
            pass
        return payload

    except subprocess.TimeoutExpired:
        try:
            logger.warning("STAR索引构建超时（120分钟）")
        except Exception:
            pass
        return {"success": False, "error": "STAR索引构建超时"}
    except Exception as e:
        try:
            logger.error(f"构建STAR索引异常：{e}")
        except Exception:
            pass
        return {"success": False, "error": f"构建STAR索引失败: {str(e)}"}


@tool
def build_hisat2_index(
    genome_id: str,
    p: Optional[int] = None,
    force_rebuild: bool = False,
) -> Dict[str, Any]:
    """构建 HISAT2 索引（等价 build_star_index）"""
    try:
        tools_config = get_tools_config()

        # 1) 从基因组配置获取基因组信息
        genomes_config_path = tools_config.genomes_config_path
        if not genomes_config_path.exists():
            return {"success": False, "error": f"基因组配置文件不存在: {genomes_config_path}"}
        
        with open(genomes_config_path, 'r', encoding='utf-8') as f:
            genome_configs = json.load(f)
        
        if genome_id not in genome_configs:
            return {"success": False, "error": f"基因组ID '{genome_id}' 未找到，可用ID: {list(genome_configs.keys())}"}
        
        genome_config = genome_configs[genome_id]
        fasta_path = genome_config.get("fasta_path")
        gtf_path = genome_config.get("gtf_path", "")  # GTF可选
        
        if not fasta_path:
            return {"success": False, "error": f"基因组 '{genome_id}' 配置缺少fasta_path"}

        # 2) 直接使用fasta路径，不做复杂转换
        fasta_file = Path(fasta_path)
        gtf_file = Path(gtf_path) if gtf_path else None
        
        # 检查必需文件是否存在
        if not fasta_file.exists():
            return {"success": False, "error": f"FASTA文件不存在: {fasta_file}，请先下载基因组文件"}
        if gtf_file and not gtf_file.exists():
            return {"success": False, "error": f"GTF文件不存在: {gtf_file}，请先下载基因组文件"}

        # 3) 直接基于fasta父目录确定索引目录
        index_dir = fasta_file.parent / "hisat2_index"
        index_dir.mkdir(parents=True, exist_ok=True)

        # 4) 检查是否需要重建
        index_files = list(index_dir.glob("genome.*.ht2"))
        if index_files and not force_rebuild:
            return {
                "success": True,
                "hisat2_index_dir": str(index_dir),
                "status": "已存在",
                "index_files": [str(f) for f in index_files],
                "message": f"HISAT2索引已存在，跳过构建: {index_dir}"
            }

        # 5) 执行索引构建
        work_dir_name = f"hisat2_index_{genome_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        work_dir = tools_config.settings.data_dir / "work" / work_dir_name
        work_dir.mkdir(parents=True, exist_ok=True)

        nf_params = {
            "genome_fasta": str(fasta_file),
            "genome_gtf": str(gtf_file) if gtf_file else "",
            "hisat2_index_dir": str(index_dir),
            "index_basename": "genome",
            "p": p or 4,
        }

        params_file = work_dir / "build_hisat2_index_params.json"
        with open(params_file, "w", encoding="utf-8") as f:
            json.dump(nf_params, f, indent=2, ensure_ascii=False)

        # 6) 定位并执行 Nextflow
        nextflow_script = tools_config.settings.nextflow_scripts_dir / "build_hisat2_index.nf"
        if not nextflow_script.exists():
            return {
                "success": False,
                "error": "未找到 build_hisat2_index.nf",
                "searched": [str(tools_config.settings.nextflow_scripts_dir / "build_hisat2_index.nf")]
            }

        logger.info(f"构建HISAT2索引 - 基因组: {genome_id}")
        logger.info(f"FASTA: {fasta_file}")
        logger.info(f"索引目录: {index_dir}")
        
        cmd = [
            "nextflow", "run", str(nextflow_script),
            "-params-file", str(params_file),
            "-work-dir", str(work_dir),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400, cwd=tools_config.settings.project_root)

        # 7) 检查结果
        final_index_files = list(index_dir.glob("genome.*.ht2"))
        success = result.returncode == 0 and len(final_index_files) > 0

        payload = {
            "success": success,
            "hisat2_index_dir": str(index_dir),
            "genome_id": genome_id,
            "index_files": [str(f) for f in final_index_files],
            "work_dir": str(work_dir),
            "params_file": str(params_file),
        }

        if success:
            payload["status"] = "构建成功"
            payload["message"] = f"HISAT2索引构建完成: {index_dir}"
        else:
            payload["status"] = "构建失败"
            payload["error"] = f"HISAT2索引构建失败，returncode: {result.returncode}"

        if get_tools_config().settings.debug_mode:
            payload.update({"stderr": result.stderr, "cmd": " ".join(cmd)})

        return payload

    except Exception as e:
        return {"success": False, "error": f"构建HISAT2索引失败: {str(e)}"}
