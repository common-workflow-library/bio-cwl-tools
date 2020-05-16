#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inputSAM
  type: File
- id: sampleName
  type: string
- id: samtoolsPath
  type: string

outputs:
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).bam

baseCommand: []
arguments:
- valueFrom: samtools view -bS $(inputs.inputSAM.path) > $(inputs.sampleName).bam
  shellQuote: false
id: samtoolsView
