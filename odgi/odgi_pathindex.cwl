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
        specs: [ https://identifiers.org/biotools/odgi ]
  DockerRequirement:
    dockerPull:  quay.io/biocontainers/odgi:0.4.1--py38h8e3bb3f_1

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
