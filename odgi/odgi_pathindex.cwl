#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi pathindex
doc: create a path index for a given graph

hints:
  SoftwareRequirement:
    packages:
      odgi:
        version: [ "0.4.1" ]
  DockerRequirement:
    #dockerPull: quay.io/biocontainers/odgi:0.3--py37h8b12597_0
    dockerImageId: odgi:latest
    dockerFile: |
      FROM python:slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install -y git nodejs npm bash cmake make g++ time 
      RUN git clone --recursive https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) && cd build && make install
      FROM python:slim
      WORKDIR /usr/local/bin/
      COPY --from=0 /usr/local/bin/odgi .

inputs:
  sparse_graph_index:
    type: File
    inputBinding:
      prefix: --idx=
      separate: false

arguments:
  - prefix: --out=
    valueFrom: $(inputs.sparse_graph_index.nameroot).og.xp
    separate: false

baseCommand: [ odgi, pathindex ]

outputs:
  indexed_paths:
    type: File
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).og.xp
