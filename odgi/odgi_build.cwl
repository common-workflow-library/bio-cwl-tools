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
      FROM python:slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install -y git nodejs npm bash cmake make g++ time 
      RUN git clone --recursive https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) && cd build && make install
      FROM python:slim
      WORKDIR /usr/local/bin/
      COPY --from=0 /usr/local/bin/odgi .

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
