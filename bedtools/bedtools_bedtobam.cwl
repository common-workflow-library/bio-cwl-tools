#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Convert reads in BED format to BAM.

requirements:
  InlineJavascriptRequirement: {}
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
        version: [ "2.25.0" ]
  
baseCommand: ["bedtools", "bedtobam"]
stdout: $(inputs.bed.nameroot).bam      

inputs:
  bed:
    type: File
    inputBinding:
      position: 1
      prefix: "-i"
  reference_info:
    type: File
    inputBinding:
      position: 2
      prefix: "-g"
    
outputs:
  bam:
    type: stdout
