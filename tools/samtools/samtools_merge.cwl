#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Merge multiple BAM files.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "merge"]

inputs:
  - id: output_name
    doc: name of merged bam file
    type: string
    inputBinding:
      position: 1
  - id: bams
    doc: bam files to be merged
    type:
      type: array
      items: File
    inputBinding:
      position: 2

outputs:
  - id: bam_merged
    type: File
    outputBinding:
      glob: $(inputs.output_name)
    
