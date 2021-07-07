#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: 
    - $(inputs.ReferenceGenome)
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0

inputs:

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
    gatk IndexFeatureFile -F $(inputs.mutantVCF.path) | gatk VariantFiltration -R $(inputs.ReferenceGenome.path) -V $(inputs.mutantVCF.path) -window 35 -cluster 3 -O $(inputs.sampleName).mutantfilter.vcf
  shellQuote: false

