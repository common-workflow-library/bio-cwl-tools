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
- id: ReferenceGenome
  type: File
- id: bcftoolsPath
  type: string

outputs:
- id: rawVCF
  type: File
  outputBinding:
    glob: $(inputs.sampleName).vcf

baseCommand: []
arguments:
- valueFrom: |-
    bcftools mpileup -f $(inputs.ReferenceGenome.path) $(inputs.inputBAM.path) | bcftools call -mv -Ob -o $(inputs.sampleName).vcf
  shellQuote: false
id: bcftools
