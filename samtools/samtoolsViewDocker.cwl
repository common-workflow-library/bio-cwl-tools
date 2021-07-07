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
- id: inputSAM
  type: File
- id: sampleName
  type: string
- id: docker
  type: string

outputs:
- id: rawBAM
  type: File
  outputBinding:
    glob: $(inputs.sampleName).bam

baseCommand: []
arguments:
- valueFrom: samtools view -bS $(inputs.inputSAM.path) > $(inputs.sampleName).bam
  shellQuote: false
id: samtoolsViewDocker
