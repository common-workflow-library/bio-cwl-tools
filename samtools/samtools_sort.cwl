#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Sort a bam file by read names.

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 15000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

baseCommand: ["samtools", "sort"]
arguments:
  - valueFrom: $(runtime.cores)
    prefix: -@
  - prefix: -m
    valueFrom: ${ return(parseInt(runtime.ram/runtime.cores-100).toString() + "M") }
    position: 1
    # specifies the allowed maximal memory usage per thread before
    # samtools start to outsource memory to temporary files

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
  
  
