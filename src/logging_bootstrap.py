"""RNA-seq 智能分析助手 · 统一日志配置模块

基于docs/debug-observability.md实现统一的日志管理：
- 支持环境变量控制的分级日志
- 文件轮转与控制台双输出
- RNA-seq特定命名空间配置
- 调试模式与生产模式切换
"""

import logging
import os
import sys
from pathlib import Path
from logging.handlers import RotatingFileHandler
from typing import Optional, Any
import json


def setup_logging(
    data_dir: Optional[Path] = None,
    log_level: Optional[str] = None,
    enable_console: bool = True,
    enable_file: bool = True
) -> None:
    """设置统一的日志配置
    
    Args:
        data_dir: 数据目录路径，默认使用Settings类的data_dir配置
        log_level: 日志级别，默认从环境变量LOG_LEVEL获取或INFO
        enable_console: 是否启用控制台输出
        enable_file: 是否启用文件输出
    """
    # 获取日志级别
    if log_level is None:
        log_level = os.getenv("LOG_LEVEL", "INFO").upper()
    level = getattr(logging, log_level, logging.INFO)
    
    # 获取数据目录 - 使用Settings类统一路径配置
    if data_dir is None:
        try:
            from .config.settings import Settings
            settings = Settings()
            data_dir = settings.data_dir
        except ImportError:
            # 回退到容器工作目录，避免循环导入
            data_dir = Path(".")
    
    # 创建日志目录
    logs_dir = data_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    # 日志格式
    fmt = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    formatter = logging.Formatter(fmt)
    
    # 清理根日志处理器
    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()
    
    # 控制台处理器
    if enable_console:
        console = logging.StreamHandler(sys.stdout)
        console.setLevel(level)
        console.setFormatter(formatter)
        root.addHandler(console)
    
    # 文件处理器（5MB轮转，保留3个备份）
    if enable_file:
        file_handler = RotatingFileHandler(
            str(logs_dir / "app.log"),
            maxBytes=5 * 1024 * 1024,  # 5MB
            backupCount=3
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)
    
    # 配置RNA-seq特定命名空间
    rna_namespaces = [
        "rna.app",
        "rna.route", 
        "rna.tools",
        "rna.star",
        "rna.hisat2",
        "rna.featurecounts",
        "rna.analysis",
        "rna.graph",
        "rna.nodes"
    ]
    
    for namespace in rna_namespaces:
        logger = logging.getLogger(namespace)
        logger.setLevel(level)
        # 继承根处理器，无需单独添加处理器
    
    # 静默第三方库的详细日志（除非DEBUG模式）
    if level > logging.DEBUG:
        logging.getLogger("httpx").setLevel(logging.WARNING)
        logging.getLogger("urllib3").setLevel(logging.WARNING)
        logging.getLogger("requests").setLevel(logging.WARNING)
        logging.getLogger("langchain").setLevel(logging.INFO)
        logging.getLogger("langgraph").setLevel(logging.INFO)


def get_logger(name: str) -> logging.Logger:
    """获取指定命名空间的日志记录器
    
    Args:
        name: 日志记录器名称，建议使用rna.*前缀
        
    Returns:
        配置好的日志记录器实例
    """
    return logging.getLogger(name)


def is_debug_enabled() -> bool:
    """检查是否启用调试模式
    
    Returns:
        True如果DEBUG环境变量为true，否则False
    """
    return os.getenv("DEBUG", "false").lower() in ("true", "1", "yes")


def _llm_logging_enabled(logger: logging.Logger) -> bool:
    """判断是否允许输出LLM payload预览（info级别输出，受开关控制）。"""
    try:
        flag = os.getenv("LLM_LOG_PAYLOAD", "true").lower() in ("true", "1", "yes")
        return bool(flag)
    except Exception:
        return False


def safe_preview(obj: Any, limit: Optional[int] = None) -> str:
    """安全生成对象的可读预览字符串。
    优先 dict()/model_dump() → JSON dumps → str()；应用长度截断；任何异常均吞掉。
    """
    limit = int(os.getenv("LLM_LOG_LIMIT", str(limit or 1200)))
    text = ""
    try:
        data = None
        # Pydantic v2
        if hasattr(obj, "model_dump") and callable(getattr(obj, "model_dump")):
            data = obj.model_dump()
        # Pydantic v1
        elif hasattr(obj, "dict") and callable(getattr(obj, "dict")):
            data = obj.dict()
        # LangGraph/Result dict-like
        elif isinstance(obj, dict):
            data = obj
        else:
            # 尝试从常见容器结构中提取
            data = obj

        if data is not None:
            try:
                text = json.dumps(data, ensure_ascii=False, indent=2)
            except Exception:
                text = str(data)
        else:
            text = str(obj)
    except Exception:
        try:
            text = str(obj)
        except Exception:
            text = "<unprintable>"

    if limit and len(text) > limit:
        return text[:limit] + "... (truncated)"
    return text


def log_llm_preview(logger: logging.Logger, scope: str, payload: Any) -> None:
    """以 info 级别输出 LLM 返回内容的限长预览。失败安全，绝不抛出异常。"""
    try:
        if not _llm_logging_enabled(logger):
            return
        preview = safe_preview(payload)
        logger.info(f"LLM输出预览({scope}): {preview}")
    except Exception:
        # 绝不影响主流程
        pass


def log_startup_info() -> None:
    """记录应用启动信息"""
    logger = get_logger("rna.app")
    logger.info("RNA-seq智能分析助手启动")
    logger.info(f"日志级别: {logging.getLevelName(logger.level)}")
    logger.info(f"调试模式: {'开启' if is_debug_enabled() else '关闭'}")
    # 获取真实的数据目录路径
    try:
        from .config.settings import Settings
        settings = Settings()
        data_dir_actual = settings.data_dir
        logger.info(f"数据目录: {data_dir_actual}")
    except ImportError:
        # 循环导入时的备用处理，不重复打印数据目录
        logger.debug("Settings模块循环导入，使用默认数据目录")


if __name__ == "__main__":
    # 测试日志配置
    setup_logging()
    log_startup_info()
    
    # 测试各个级别的日志
    logger = get_logger("rna.test")
    logger.debug("这是DEBUG级别消息")
    logger.info("这是INFO级别消息")
    logger.warning("这是WARNING级别消息")
    logger.error("这是ERROR级别消息")
