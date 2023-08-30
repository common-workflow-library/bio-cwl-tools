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
        specs: [ https://identifiers.org/biotools/odgi ]
  DockerRequirement:
    dockerPull:  quay.io/biocontainers/odgi:0.4.1--py38h8e3bb3f_1

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
