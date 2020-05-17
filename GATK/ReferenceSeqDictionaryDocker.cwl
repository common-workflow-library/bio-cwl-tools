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
    java -jar $(inputs.gatk_path) CreateSequenceDictionary -R $(inputs.ReferenceGenome.path) -O $(inputs.sampleName).fasta.dict
  shellQuote: false
id: ReferenceSeqDictionaryDocker
