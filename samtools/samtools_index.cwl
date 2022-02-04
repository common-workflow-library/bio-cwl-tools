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
    dockerPull: quay.io/biocontainers/samtools:1.14--hb421002_0

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
    format: edam:format_2572  # BAM 
    outputBinding:
      glob: $(inputs.bam_sorted.basename)
      
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
