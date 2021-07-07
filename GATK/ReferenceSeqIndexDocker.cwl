#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: CommandLineTool
cwlVersion: v1.0

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

inputs:
- id: ReferenceGenome
  type: File
- id: sampleName
  type: string

outputs:
- id: refIndex
  type: File
  outputBinding:
    glob: $(inputs.sampleName).fasta.fai

baseCommand: []
arguments:
- valueFrom: samtools faidx $(inputs.ReferenceGenome.path) > $(inputs.sampleName).fasta.fai
  shellQuote: false
id: ReferenceSeqIndexDocker
