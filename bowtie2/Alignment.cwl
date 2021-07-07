#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: leftFastq
  type: File
- id: rightFastq
  type: File
- id: ReferenceGenome
  type: File
- id: sampleName
  type: string
- id: index
  type: string

outputs:
- id: rawSAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).sam

baseCommand: []
arguments:
- valueFrom: |-
    export PATH=$PATH:/home/ngsap2/Downloads/bowtie2-2.4.1 | bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index) | bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam
  shellQuote: false
id: Alignment
