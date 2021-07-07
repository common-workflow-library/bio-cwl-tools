#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: GATK
  type: File
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
    $(inputs.GATK.path) CreateSequenceDictionary -R $(inputs.ReferenceGenome.path) -O $(inputs.sampleName).fasta.dict
  shellQuote: false
id: ReferenceSeqDictionary
