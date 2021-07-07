#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: 
    - $(inputs.Reference2bitGenome)
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0 

inputs:

- id: inputBAM
  type: File
- id: Reference2bitGenome
  type: File
- id: RefIndex
  type: File
- id: RefDict
  type: File
- id: sampleName
  type: string
- id: BAMindex
  type: File
- id: output_filename
  type: string

outputs:
- id: rawVCF
  type: File
  outputBinding:
    glob: $(inputs.output_filename)

baseCommand: []
arguments:
- valueFrom: |-
    gatk HaplotypeCallerSpark -R $(inputs.Reference2bitGenome.path) -I $(inputs.inputBAM.path) -O $(inputs.output_filename)
  shellQuote: false

