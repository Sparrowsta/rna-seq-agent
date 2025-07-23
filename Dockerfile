# 1. 使用官方Ubuntu 22.04作为基础镜像
FROM ubuntu:22.04

# 2. 设置一些环境变量，避免安装过程中的交互提示
ENV DEBIAN_FRONTEND=noninteractive

# 3. 更换apt源为清华源，然后更新包管理器并安装一些基础工具
RUN sed -i 's/archive.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list && \
    sed -i 's/security.ubuntu.com/mirrors.tuna.tsinghua.edu.cn/g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends wget ca-certificates procps python3 python3-pip && \
    apt-get install -y openjdk-21-jdk-headless && \
    apt-get clean && \
    pip3 install nextflow && \
    rm -rf /var/lib/apt/lists/*

# 4. 下载并安装Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# 5. 将conda添加到PATH，并配置清华镜像源
ENV PATH=$CONDA_DIR/bin:$PATH
RUN ln -s $CONDA_DIR/bin/conda /usr/local/bin/conda && \
    ln -s $CONDA_DIR/bin/activate /usr/local/bin/activate && \
    ln -s $CONDA_DIR/bin/deactivate /usr/local/bin/deactivate && \
    # to ensure the solver can find compatible packages efficiently.
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --set show_channel_urls yes

# 6. 接受服务条款，然后为每个工具或工具组创建独立的Conda环境
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda install -y -c conda-forge mamba && \
    mamba create -y -n sra_env -c bioconda sra-tools && \
    mamba create -y -n qc_env -c bioconda fastp=1.0.1 && \
    mamba create -y -n align_env -c bioconda samtools=1.22.1 star=2.7.11b && \
    mamba create -y -n quant_env -c bioconda subread=2.1.1 && \
    # Create a dedicated environment for our Python scripts and their dependencies
    mamba create -y -n ngs_env python=3.12 && \
    # Install Python packages using mamba from the conda-forge channel for better dependency management.
    mamba install -y -n ngs_env -c conda-forge pandas python-dotenv langchain langchain-openai tabulate requests && \
    mamba clean -y -a

# 7. 创建工作目录并将项目文件复制到镜像中
WORKDIR /app
COPY . .

# 8. 设置入口点，使容器可以直接运行启动脚本

ENTRYPOINT ["/bin/sh", "-c", "exec python3 launch.py"]
