#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: odgi viz
doc: variation graph visualizations

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: 1
  SoftwareRequirement:
    packages:
      odgi:
        version: [ "0.4.1" ]
        specs: [ https://identifiers.org/biotools/odgi ]
  DockerRequirement:
    dockerPull:  quay.io/biocontainers/odgi:0.4.1--py38h8e3bb3f_1

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
  - --out=$(inputs.sparse_graph_index.nameroot).png

baseCommand: [ odgi, viz ]

outputs:
  graph_image:
    type: File
    format: iana:image/png
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).png

$namespaces:
  iana: https://www.iana.org/assignments/media-types/
