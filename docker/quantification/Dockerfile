FROM continuumio/miniconda3

WORKDIR /opt


RUN apt-get update \
    && apt-get -y install curl \
    && apt-get -y install wget \
    && apt-get -y install git\
    && apt-get -y install \
    python3-pip \
 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

RUN curl https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/v1_index.tar.gz | tar -xzf -
RUN curl https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/v2_index.tar.gz | tar -xzf -
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visiumv1.tx2gene.tsv
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visiumv2.tx2gene.tsv
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v1.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v2.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v3.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v4.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v5.txt

COPY docker/quantification/requirements.txt /opt/requirements.txt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN conda install -c bioconda samtools

RUN mkdir /out

COPY bin/common /opt
COPY bin/quantification /opt
