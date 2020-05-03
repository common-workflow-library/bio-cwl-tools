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
      FROM debian:stable-slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install -y git bash cmake make g++ python-dev
      RUN git clone --recursive https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -H. -Bbuild && cmake --build build -- -j $(nproc) && cd build && make install
      FROM debian:stable-slim
      RUN apt-get update && apt-get install -y libgomp1 && rm -rf /var/lib/apt/lists/*
      COPY --from=0 /usr/local/bin/odgi /usr/local/bin/

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
  - --fasta=$(inputs.sparse_graph_index.nameroot).og.fasta

stdout: $(inputs.sparse_graph_index.nameroot).w$(inputs.bin_width).json

baseCommand: [ odgi, bin ]

outputs:
  bins:
    type: stdout
    format: iana:application/json
  pangenome_sequence:
    type: File
    format: edam:format_1929  # FASTA
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).og.fasta

$namespaces:
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
  - http://edamontology.org/EDAM_1.18.owl
