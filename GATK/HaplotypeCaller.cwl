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
    glob: $(inputs.sampleName).mutant.vcf

baseCommand: []
arguments:
- valueFrom: |-
    samtools index $(inputs.inputBAM.path) > $(inputs.sampleName).split.bam.bai                                      $(inputs.GATK.path) HaplotypeCaller -R $(inputs.ReferenceGenome.path) -I $(inputs.inputBAM.path) -O $(inputs.sampleName).mutant.vcf
  shellQuote: false
id: HaplotypeCaller
