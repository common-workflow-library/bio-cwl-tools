#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: 
    - $(inputs.inputBAM)
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0 

inputs:

- id: inputBAM
  type: File
- id: sampleName
  type: string

outputs:
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).bam

baseCommand: []
arguments:
- valueFrom: |-
    gatk SortSamSpark -I $(inputs.inputBAM.path) -O $(inputs.sampleName).bam
  shellQuote: false
id: SortSamSpark
