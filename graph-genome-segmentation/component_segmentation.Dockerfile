FROM debian:stable-slim
RUN apt-get update && apt-get install -y git python3-pip && rm -rf /var/lib/apt/lists/*
WORKDIR /src
RUN git clone --depth 1 -b optimizations https://github.com/graph-genome/component_segmentation && \
  rm -Rf component_segmentation/data
WORKDIR /src/component_segmentation
RUN pip3 install --no-cache-dir -r requirements.txt
RUN chmod a+x segmentation.py
ENV PATH="/src/component_segmentation:${PATH}"
