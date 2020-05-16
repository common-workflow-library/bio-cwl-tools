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
- id: inputBAM
  type: File
- id: ReferenceGenome
  type: File
- id: sampleName
  type: string
- id: docker
  type: string

outputs:
- id: rawMpileup
  type: File
  outputBinding:
    glob: $(inputs.sampleName).mpileup

baseCommand: []
arguments:
- valueFrom: |-
    samtools index $(inputs.inputBAM.path)samtools mpileup -B -f $(inputs.ReferenceGenome.path) $(inputs.inputBAM.path) > $(inputs.sampleName).mpileup
  shellQuote: false
id: samtoolsMpileupDocker
