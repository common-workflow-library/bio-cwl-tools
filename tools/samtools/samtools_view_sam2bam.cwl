#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Convert SAM to BAM.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "view"]
arguments:
  - valueFrom: -h
    position: 1
    # include the headers
  - valueFrom: -b
    position: 1
    # output in bam format
stdout: $(inputs.sam.nameroot).bam

inputs:
  sam:
    doc: reads to be checked in sam format
    type: File
    inputBinding:
      position: 2

outputs:
  bam:
    type: stdout
  
  
