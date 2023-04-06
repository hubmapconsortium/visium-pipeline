
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    openjdk-17-jre-headless \
    python3 \
    unzip \
 && rm -rf /var/cache/apt/*

ADD https://downloads.openmicroscopy.org/bio-formats/6.10.0/artifacts/bftools.zip /opt/bftools.zip

WORKDIR /opt
RUN unzip bftools.zip \
 && rm bftools.zip

COPY bin/bftools_wrapper.py /opt