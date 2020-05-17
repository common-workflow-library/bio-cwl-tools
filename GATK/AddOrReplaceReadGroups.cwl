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
- id: inputSAM
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
    $(inputs.GATK.path) AddOrReplaceReadGroups -I $(inputs.inputSAM.path) -O $(inputs.sampleName).bam -RGID 1 -RGLB 445_LIB -RGPL illumina -RGSM RNA -RGPU illumina
  shellQuote: false
id: AddOrReplaceReadGroups
