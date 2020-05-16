#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inputBAM
  type: File
- id: sampleName
  type: string
- id: samtoolsPath
  type: string

outputs:
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).bam.bai

baseCommand: []
arguments:
- valueFrom: samtools index $(inputs.inputBAM.path) > $(inputs.sampleName).bam.bai
  shellQuote: false
id: samtoolsIndex
