#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi build
doc: construct a dynamic succinct variation graph

requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.graph)
        writable: true

hints:
  SoftwareRequirement:
    packages:
      odgi:
        version: [ "0.4.1" ]
  DockerRequirement:
    #dockerPull: quay.io/biocontainers/odgi:0.3--py37h8b12597_0
    dockerImageId: odgi:latest
    dockerFile: |
      FROM debian:stable-slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install -y git bash cmake make g++ python-dev && rm -rf /var/lib/apt/lists/*
      RUN git clone --recursive https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) && cd build && make install
      FROM debian:stable-slim
      RUN apt-get update && apt-get install -y libgomp1 && rm -rf /var/lib/apt/lists/*
      COPY --from=0 /usr/local/bin/odgi /usr/local/bin/

inputs:
  graph:
    type: File
    inputBinding:
      prefix: --gfa=
      separate: false
      valueFrom: $(self.basename)
    #format: GFA

arguments:
  - --progress
  - prefix: --out=
    valueFrom: $(inputs.graph.nameroot).og
    separate: false

baseCommand: [ odgi, build ]

outputs:
  sparse_graph_index:
    type: File
    outputBinding:
      glob: $(inputs.graph.nameroot).og  
