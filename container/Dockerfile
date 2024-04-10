FROM ubuntu:22.04

WORKDIR /scripts

COPY . .

RUN apt-get update \
    && apt-get install -yq --no-install-recommends \
	build-essential zlib1g-dev libbz2-dev liblzma-dev unzip wget python3 python3-pip \
    && pip3 install pandas

RUN make install

# adding the path to the tools managed by mamba
ENV PATH="/scripts/tools/bin:$PATH"

CMD ["bash","/scripts/liftover_hg19_to_hg38.sh"]