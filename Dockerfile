FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        g++ \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install anaconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair3 python=3.9.0 -y

ENV PATH /opt/conda/envs/clair3/bin:$PATH
ENV CONDA_DEFAULT_ENV clair3

RUN /bin/bash -c "source activate clair3" && \
    conda install -c conda-forge pypy3.6 -y && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install mpmath==1.2.1 && \
    conda install -c conda-forge tensorflow==2.8.0 -y && \
    conda install -c conda-forge pytables -y && \
    conda install -c anaconda pigz cffi==1.14.4 -y && \
    conda install -c conda-forge parallel=20191122 zstd -y && \
    conda install -c conda-forge -c bioconda samtools=1.15.1 -y && \
    conda install -c conda-forge -c bioconda whatshap=1.7 -y && \
    conda install -c conda-forge xz zlib bzip2 automake curl -y && \
    pip install tensorflow-addons && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip && \
    echo "source activate clair3" > ~/.bashrc

COPY . .

RUN wget http://www.bio8.cs.hku.hk/clair3_trio/clair3_models/clair3_models.tar.gz -P /opt/models/clair3_models && \
    tar -zxvf /opt/models/clair3_models/clair3_models.tar.gz -C /opt/models/clair3_models && \
    rm /opt/models/clair3_models/clair3_models.tar.gz && \
	wget http://www.bio8.cs.hku.hk/clair3_trio/clair3_trio_models/clair3_trio_models.tar.gz -P /opt/models/clair3_trio_models && \
    tar -zxvf /opt/models/clair3_trio_models/clair3_trio_models.tar.gz -C /opt/models/clair3_trio_models && \
    rm /opt/models/clair3_trio_models/clair3_trio_models.tar.gz && \
    cd /opt/bin && \
    PREFIX=/opt/conda/envs/clair3 PYTHON=/opt/conda/envs/clair3/bin/python && \
    rm -rf /opt/bin/samtools-* /opt/bin/longphase-*
