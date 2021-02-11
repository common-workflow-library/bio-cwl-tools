FROM debian:stable-slim
RUN apt-get update && apt-get install -y git python3-pip && rm -rf /var/lib/apt/lists/*
WORKDIR /src
RUN git clone --depth 10 https://github.com/graph-genome/component_segmentation && \
  rm -Rf component_segmentation/data
WORKDIR /src/component_segmentation
RUN git checkout 22f16d5669e956eada41e1c40ca0ecb7ea9823aa
RUN pip3 install --no-cache-dir -r requirements.txt .
