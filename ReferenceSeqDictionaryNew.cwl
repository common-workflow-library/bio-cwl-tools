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

- id: ReferenceGenome
  type: File
- id: sampleName
  type: string

outputs:
- id: refDict
  type: File
  outputBinding:
    glob: $(inputs.sampleName).fasta.dict

baseCommand: []
arguments:
- valueFrom: |-
    gatk CreateSequenceDictionary -R $(inputs.ReferenceGenome.path) -O $(inputs.sampleName).fasta.dict
  shellQuote: false
id: ReferenceSeqDictionary
