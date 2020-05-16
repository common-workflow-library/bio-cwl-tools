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
- id: inputMpileup
  type: File
- id: sampleName
  type: string
- id: docker
  type: string

outputs:
- id: snpVCF
  type: File
  outputBinding:
    glob: $(inputs.sampleName).vcf

baseCommand: []
arguments:
- valueFrom: varscan mpileup2snp $(inputs.inputMpileup.path) > $(inputs.sampleName).vcf
  shellQuote: false
id: mpileup2snpDocker
