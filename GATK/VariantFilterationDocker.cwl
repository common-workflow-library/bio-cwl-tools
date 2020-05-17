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
- id: docker
  type: string
- id: gatk_path
  type: string
- id: mutantVCF
  type: File
- id: ReferenceGenome
  type: File
- id: RefIndex
  type: File
- id: RefDict
  type: File
- id: sampleName
  type: string

outputs:
- id: rawVCF
  type: File
  outputBinding:
    glob: $(inputs.sampleName).mutantfilter.vcf

baseCommand: []
arguments:
- valueFrom: |-
    java -jar $(inputs.gatk_path) IndexFeatureFile -F $(inputs.mutantVCF.path)      java -jar $(inputs.gatk_path) VariantFiltration -R $(inputs.ReferenceGenome.path) -V $(inputs.mutantVCF.path) -window 35 -cluster 3 -O $(inputs.sampleName).mutantfilter.vcf
  shellQuote: false
id: VariantFilterationDocker
