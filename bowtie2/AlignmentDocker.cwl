#!/usr/bin/env cwl-runner

# This tool description was generated automatically by wdl2cwl ver. 0.2



class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: docker

inputs:
  ReferenceGenome:
    id: ReferenceGenome
    type: File
  docker:
    id: docker
    type: string
  index:
    id: index
    type: string
  leftFastq:
    id: leftFastq
    type: File
  rightFastq:
    id: rightFastq
    type: File
  sampleName:
    id: sampleName
    type: string

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.sampleName).sam

baseCommand: []
arguments:
  valueFrom: |-
    bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index) | bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam
  shellQuote: false
