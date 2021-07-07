#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: inputMpileup
  type: File
- id: sampleName
  type: string
- id: varscanPath
  type: string

outputs:
- id: snpVCF
  type: File
  outputBinding:
    glob: $(inputs.sampleName).vcf

baseCommand: []
arguments:
- valueFrom: |-
    java -jar $(inputs.varscanPath) mpileup2snp $(inputs.inputMpileup.path) > $(inputs.sampleName).vcf
  shellQuote: false
id: mpileup2snp
