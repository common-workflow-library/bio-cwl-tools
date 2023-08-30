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
        specs: [ https://identifiers.org/biotools/odgi ]
  DockerRequirement:
    dockerPull:  quay.io/biocontainers/odgi:0.4.1--py38h8e3bb3f_1

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
  - https://edamontology.org/EDAM_1.18.owl
