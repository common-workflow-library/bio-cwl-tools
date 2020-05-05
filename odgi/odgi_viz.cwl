#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi viz
doc: variation graph visualizations

requirements:
  InlineJavascriptRequirement: {}
hints:
  SoftwareRequirement:
    packages:
      odgi:
        version: [ "0.4.1" ]
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: 1
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
  sparse_graph_index: File
  width:
    type: int?
    doc: width in pixels of output image
    inputBinding:
      prefix: --width=
      separate: false
  height:
    type: int?
    doc: height in pixels of output image
    inputBinding:
      prefix: --height=
      separate: false
  path_per_row:
    type: boolean?
    doc: display a single path per row rather than packing them
    inputBinding:
      prefix: --path-per-row
  path_height:
    type: int?
    doc: path display height
    inputBinding:
      prefix: --path-height=
      separate: false

arguments:
  - --idx=$(inputs.sparse_graph_index.path)
  - --threads=$(runtime.cores)
  - --out=$(inputs.sparse_graph.index.nameroot).png

baseCommand: [ odgi, viz ]

outputs:
  graph_image:
    type: File
    format: iana:image/png
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).png

$namespaces:
  iana: https://www.iana.org/assignments/media-types/
