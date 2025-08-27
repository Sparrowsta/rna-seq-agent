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

# 根据用户要求，设置所有Nextflow相关路径的环境变量
# 设置主目录以解决用户上下文和权限问题

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
        python3-pip \
        openjdk-21-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# 5. 安装Nextflow
RUN pip3 install --no-cache-dir nextflow

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

# 8. 安装Python依赖
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt


# 9. 设置根目录为工作目录，数据通过volume映射

# 10. 设置工作目录为/data，与数据映射一致
WORKDIR /data

# 11. 复制应用文件到根目录（应用与数据分离）
COPY main.nf /
COPY main.py /
COPY src/ /src/

# 12. 复制静态配置文件到镜像
COPY config/genomes.json /config/genomes.json

# data目录和.env文件将通过volume映射和环境变量传入

# 文件权限将由运行时的 --user 参数控制

# 13. 设置环境变量
ENV HOME=/data \
    NXF_HOME=/data/.nextflow \
    NXF_WORK=/data/work \
    NXF_TEMP=/data/tmp \
    PYTHONPATH=/src

# 14. 健康检查
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python3 -c "import sys; sys.exit(0)"

# 15. 启动命令
CMD ["python3", "-u", "/main.py"]
