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
    glob: $(inputs.sampleName).sorted.bam

baseCommand: []
arguments:
- valueFrom: samtools sort -l 0 -o $(inputs.sampleName).sorted.bam $(inputs.inputBAM.path)
  shellQuote: false
id: samtoolsSort
