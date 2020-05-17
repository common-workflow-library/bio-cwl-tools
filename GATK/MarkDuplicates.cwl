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
- id: inputBAM
  type: File
- id: sampleName
  type: string

outputs:
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).markdup.bam

baseCommand: []
arguments:
- valueFrom: |-
    $(inputs.GATK.path) MarkDuplicates -I $(inputs.inputBAM.path) -O $(inputs.sampleName).markdup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -M output.metrics
  shellQuote: false
id: MarkDuplicates
