#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/varscan:2.3.7--2

inputs:
- id: inputMpileup
  type: File
- id: sampleName
  type: string

outputs:
- id: snpVCF
  type: File
  outputBinding:
    glob: $(inputs.sampleName).snp.vcf

baseCommand: []
arguments:
- valueFrom: varscan mpileup2snp $(inputs.inputMpileup.path) > $(inputs.sampleName).snp.vcf
  shellQuote: false

