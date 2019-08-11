#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Indexing BAM.

requirements:
  InitialWorkDirRequirement:
    listing: 
      - $(inputs.bam_sorted)
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "index"]
arguments:
  - valueFrom: -b  # specifies that index is created in bai format
    position: 1

inputs:
  bam_sorted:
    doc: sorted bam input file
    type: File
    inputBinding:
      position: 2

outputs:
  bam_sorted_indexed:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.bam_sorted.basename)
      
    
