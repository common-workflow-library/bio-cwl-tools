#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "broadinstitute/gatk:4.1.3.0"
  InlineJavascriptRequirement: {}

inputs:
  # REQUIRED ARGS

  Variant:
    type: File
    inputBinding:
      prefix: "-V"

  Resource:
    type: File
    inputBinding:
      prefix: "-resource"
    secondaryFiles:
      - .idx

  Output:
    type: string
    default: "filtered.vcf"
    inputBinding:
      prefix: "-O"
      valueFrom: "filtered.vcf"

  # OPTIONAL ARGS
  Reference:
    type: File?
    inputBinding:
      prefix: "-R"
    secondaryFiles:
      - ^.dict
      - .fai

  Architecture:
    type: string?
    inputBinding:
      prefix: "-architecture"

  Tranche:
    type: int?
    inputBinding:
      prefix: "-tranche"

  InfoKey:
    type: string?
    inputBinding:
      prefix: "--info-key"

baseCommand: ["/gatk/gatk"]

arguments:
  - valueFrom: "FilterVariantTranches"
    position: -1

outputs:
  filteredVCF:
    type: File
    outputBinding:
      glob: "filtered.vcf"
