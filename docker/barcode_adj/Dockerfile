FROM pypy:3.8

## System packages
RUN apt-get update \
 && apt-get -y install --no-install-recommends \
    samtools \
 && rm -rf /var/cache/apt/* /var/lib/apt/lists/*

## Python packages
COPY docker/barcode_adj/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin/common /opt
COPY bin/analysis /opt

WORKDIR /opt

RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v1.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v2.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v3.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v4.txt
RUN wget https://s3.amazonaws.com/hubmap-pipeline-assets/salmon/visium/visium-v5.txt
