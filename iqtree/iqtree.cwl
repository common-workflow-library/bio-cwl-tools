#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  alignment:
    type: File
    label: Multiple Sequence Alignment

outputs:
  tree:
    type: File
    label: Maximum-likelihood tree
    outputBinding:
      glob: "*.treefile"

arguments:
  - -s 
  - $(inputs.alignment.path)

baseCommand: iqtree2

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/iqtree:2.2.2.7--h21ec9f0_2

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.alignment)
