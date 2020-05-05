#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi build
doc: construct a dynamic succinct variation graph

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.graph)
        writable: true

hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: $(Math.ceil((inputs.graph.size/(1024*1024*1024)+1) * 2))
  SoftwareRequirement:
    packages:
      odgi:
        version: [ "0.4.1" ]
  DockerRequirement:
    #dockerPull: quay.io/biocontainers/odgi:0.3--py37h8b12597_0
    dockerImageId: odgi:latest
    dockerFile: |
      FROM debian:bullseye-slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install --no-install-recommends -y \
        ca-certificates \
        bash \
        cmake \
        git \
        g++ \
        make \
        python-dev \
        && rm -rf /var/lib/apt/lists/*
      RUN git clone --recursive --branch v0.4.1 https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -DCMAKE_BUILD_TYPE=Release -H. -Bbuild && \
        cmake --build build --config Release -- -j $(nproc) && \
        cmake --install build/ --config Release -v --strip
      FROM debian:bullseye-slim
      RUN apt-get update && apt-get install -y libgomp1 && rm -rf /var/lib/apt/lists/*
      COPY --from=0 /usr/local/bin/odgi /usr/local/bin/

inputs:
  graph:
    type: File
    #format: GFA1 or GFA2

  sort:
    type: boolean?
    doc: apply generalized topological sort to the graph and set node ids to order
    inputBinding:
      prefix: --sort

arguments:
  - --progress
  - --gfa=$(inputs.graph.basename)
  - --out=-

stdout: $(inputs.graph.nameroot).odgi

baseCommand: [ odgi, build ]

outputs:
  sparse_graph_index: stdout
