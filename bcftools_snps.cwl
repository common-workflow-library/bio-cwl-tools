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
  dockerPull: quay.io/biocontainers/bcftools:1.7--0


inputs:

- id: inputBAM

  type: File

- id: sampleName

  type: string

- id: ReferenceGenome

  type: File


outputs:

- id: snp

  type: File

  outputBinding:

    glob: $(inputs.sampleName).snp.vcf



baseCommand: []

arguments:

- valueFrom: |-

    bcftools mpileup -f $(inputs.ReferenceGenome.path) $(inputs.inputBAM.path) | bcftools call --skip-variants indels -mv -Ov -o $(inputs.sampleName).snp.vcf

  shellQuote: false


