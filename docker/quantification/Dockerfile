FROM ubuntu:20.04

WORKDIR /opt


RUN apt-get update \
    && apt-get -y install curl \
    && apt-get -y install wget \
    && apt-get -y install git\
    && apt-get -y install \
    python3-pip \
    && apt-get -y install libncurses5 \

 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.fai
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.bwt
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.amb
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.ann
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.sa
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_Mouse_v1.0.1.fasta.pac
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.fai
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.bwt
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.amb
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.ann
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.sa
RUN wget https://hubmap-pipeline-assets.s3.amazonaws.com/salmon/rna-probe-sets/probe-sets/Chromium_v1.0.1.fasta.pac

COPY docker/quantification/requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN conda install -c bioconda samtools
RUN conda install -c bioconda bwa


RUN mkdir /out

COPY bin/common /opt
COPY bin/quantification /opt
