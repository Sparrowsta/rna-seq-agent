# 1. 使用官方Ubuntu 22.04作为基础镜像
FROM docker.m.daocloud.io/ubuntu:22.04

# 2. 设置一些环境变量，避免安装过程中的交互提示
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# 3. 更换apt源为清华源，然后更新包管理器并安装一些基础工具
RUN sed -i 's/archive.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list && \
    sed -i 's/security.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends wget ca-certificates procps python3 python3-pip curl uvicorn jq && \
    apt-get install -y openjdk-21-jdk-headless && \
    apt-get clean && \
    pip3 install nextflow && \
    rm -rf /var/lib/apt/lists/*

# 4. 下载并安装Miniconda (使用curl和清华镜像源加速)
ENV CONDA_DIR=/opt/conda
RUN curl -L -o ~/miniconda.sh https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# 5. 将conda添加到PATH，并配置清华镜像源
ENV PATH=$CONDA_DIR/bin:$PATH
RUN ln -s $CONDA_DIR/bin/conda /usr/local/bin/conda && \
    ln -s $CONDA_DIR/bin/activate /usr/local/bin/activate && \
    ln -s $CONDA_DIR/bin/deactivate /usr/local/bin/deactivate && \
    echo "channels:" > /root/.condarc && \
    echo "  - defaults" >> /root/.condarc && \
    echo "show_channel_urls: true" >> /root/.condarc && \
    echo "default_channels:" >> /root/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main" >> /root/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r" >> /root/.condarc && \
    echo "  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2" >> /root/.condarc && \
    echo "custom_channels:" >> /root/.condarc && \
    echo "  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /root/.condarc && \
    echo "  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /root/.condarc && \
    echo "  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud" >> /root/.condarc

# 6. 为每个工具或工具组创建独立的Conda环境
RUN conda install -y -c conda-forge mamba


RUN mamba create -y -n sra_env -c conda-forge -c bioconda sra-tools=3.2.0 aspera-cli=4.20.0
RUN mamba create -y -n qc_env -c conda-forge -c bioconda fastp=1.0.1
RUN mamba create -y -n align_env -c conda-forge -c bioconda samtools=1.22.1 star=2.7.11b
RUN mamba create -y -n quant_env -c conda-forge -c bioconda subread=2.1.1


# 直接在全局 Python 环境中安装依赖
RUN pip install --no-cache-dir pandas python-dotenv 'langchain>=0.2.0' 'langchain-openai>=0.1.0' 'fastapi>=0.110.0' 'pydantic>=2.0.0' sse-starlette tabulate requests mcp langgraph langchain-google-genai 


# Create a dedicated environment for differential expression analysis with R
RUN mamba create -y -n de_env -c conda-forge -c r r-base r-essentials r-argparse

# Set up CRAN and Bioconductor mirrors for faster R package installation
# This creates a global .Rprofile that R will load on startup.
RUN echo "options(repos = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'), BioC_mirror = 'https://mirrors.tuna.tsinghua.edu.cn/bioconductor')" > /root/.Rprofile

# Install BiocManager and required Bioconductor packages
# The mirrors are now automatically picked up from the .Rprofile file.
RUN conda run --no-capture-output -n de_env R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('DESeq2', 'EnhancedVolcano', 'pheatmap'), update=FALSE, ask=FALSE)"


RUN mamba clean -y -a

# 7. 创建工作目录并精确复制项目文件
WORKDIR /app

# 7.1. 创建数据目录
# 我们使用 RUN mkdir -p 来创建 data 目录，因为它当前不存在于项目中。
# 这避免了在本地创建一个空目录再复制进来的需要，保持了构建上下文的整洁。
RUN mkdir -p /app/data

# 7.2. 复制配置文件和 Nextflow 脚本 (作为不常变动的文件)
COPY main.nf ./

COPY config/genomes.json ./config/
COPY config/nextflow.config ./config/
# 7.3. 复制应用入口文件
COPY main.py ./

# 7.4. 复制核心应用代码 (最常变动，放在最后)
COPY agent/ ./agent/

# 8. 定义容器启动时要执行的默认命令
# 直接使用全局的 python3 运行服务器脚本，彻底绕开 conda run
CMD ["python3", "-u", "main.py"]
