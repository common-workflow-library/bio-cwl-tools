#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Clips regions in a bed file that are exceeding chromosome boundaries.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
  DockerRequirement:
    dockerPull: biocontainers/bedtools:2.25.0
  SoftwareRequirement:
    packages:
      bedtools:
        specs: [ "http://identifiers.org/biotools/bedtools" ]
        version: [ "2.26.0" ]

baseCommand: ["bedtools", "slop"]
arguments:
  - valueFrom: "0"
    prefix: -b
    position: 1
stdout: $(inputs.bed.nameroot)_clipped.bed
  
inputs:
  bed:
    type: File
    inputBinding:
      prefix: "-i"
      position: 2
  reference_info:
    type: File
    inputBinding:
      prefix: "-g"
      position: 3

outputs:
  bed_clipped:
    type: stdout
    
