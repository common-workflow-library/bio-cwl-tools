#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Sort a bam file by read names.

requirements:
  InlineJavascriptRequirement: {}
hints:
  
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "sort"]

inputs:
  bam_unsorted:
    doc: aligned reads to be checked in sam or bam format
    type: File
    inputBinding:
      position: 2
  by_name:
    doc: If true, will sort by name, otherwise will sort by genomic position
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -n

stdout: $(inputs.bam_unsorted.basename)

outputs:
  bam_sorted:
    type: stdout
  
  
