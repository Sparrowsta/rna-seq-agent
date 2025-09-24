# 1. 使用官方micromamba镜像作为基础镜像
FROM mambaorg/micromamba:1.5-jammy

# 2. 切换到root用户进行系统配置
USER root

# 3. 设置环境变量，优化Python和避免交互提示
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Asia/Shanghai \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# 4.1 使用模板写入清华 APT 源（根据发行版代号自动适配）
# 基于 /etc/os-release 的 VERSION_CODENAME（本镜像为 jammy）；保留 security 官方源
RUN cat > /etc/apt/sources.list <<EOF
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse
deb http://security.ubuntu.com/ubuntu/ jammy-security main restricted universe multiverse
EOF

# 4. 安装必要的系统工具和Java
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        ca-certificates \
        procps \
        curl \
        jq \
        gosu \
        git \
        python3 \
        python3-venv \
        python3-pip \
        openjdk-21-jre-headless \
        tzdata && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 5. 安装uv（现代Python包管理器）
RUN curl -LsSf https://astral.sh/uv/install.sh | sh \
    && install -m 0755 /root/.local/bin/uv /usr/local/bin/uv
ENV PATH="/root/.local/bin:$PATH"

# 6. 配置mamba镜像源
RUN echo "channels:" > /opt/conda/.condarc && \
    echo "  - defaults" >> /opt/conda/.condarc && \
    echo "show_channel_urls: true" >> /opt/conda/.condarc && \
    echo "default_channels:" >> /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main" >> /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r" >> /opt/conda/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2" >> /opt/conda/.condarc && \
    echo "custom_channels:" >> /opt/conda/.condarc && \
    echo "  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /opt/conda/.condarc && \
    echo "  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /opt/conda/.condarc && \
    echo "  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /opt/conda/.condarc

# 7. 创建生物信息学conda环境
RUN micromamba create -y -n sra_env -c conda-forge -c bioconda sra-tools=3.2.0 aspera-cli=4.20.0 && \
    micromamba create -y -n qc_env -c conda-forge -c bioconda fastp=1.0.1 && \
    micromamba create -y -n align_env -c conda-forge -c bioconda samtools=1.22.1 star=2.7.11b hisat2=2.2.1 && \
    micromamba create -y -n quant_env -c conda-forge -c bioconda subread=2.1.1 && \
    micromamba clean -y -a

# 8. 设置工作目录为/data，与数据映射一致
WORKDIR /data

# 9. 复制项目文件
COPY pyproject.toml uv.lock ./
COPY main.py /
COPY src/ /src/

# 10. 使用uv安装Python依赖（基于pyproject.toml）
RUN uv pip install --system -r pyproject.toml

# 11. 设置环境变量
ENV HOME=/data \
    NXF_HOME=/data/.nextflow \
    PYTHONPATH=/src \
    NXF_OFFLINE=1 \
    PATH="/root/.local/bin:/root/.cargo/bin:$PATH"

# 12. 健康检查
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD uv run python -c "import sys; sys.exit(0)"

# 13. 启动命令
CMD ["uv", "run", "python", "-u", "/main.py"]
