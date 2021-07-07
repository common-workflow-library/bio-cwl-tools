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
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).split.bam

baseCommand: []
arguments:
- valueFrom: |-
    java -jar $(inputs.gatk_path) SplitNCigarReads -R $(inputs.ReferenceGenome.path) -I $(inputs.inputBAM.path) -O $(inputs.sampleName).split.bam
  shellQuote: false
id: SplitNCigarReadsDocker
