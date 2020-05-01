#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi bin
doc: binning of path information in the graph

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

  bin_width:
    type: int?
    doc: width of each bin in basepairs along the graph vector
    inputBinding:
      prefix: --bin-width=
      separate: false

arguments:
  - --json

stdout: $(inputs.sparse_graph_index.nameroot).w$(inputs.bin_width).json

baseCommand: [ odgi, bin ]

outputs:
  bins: stdout
