#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
  StepInputExpressionRequirement: {}

inputs:
  InputFile:
    type: File
    format: edam:format_1929  # FASTA
    inputBinding:
      position: 200

#Optional arguments

  algoType:
    type: 
      - "null"
      - type: enum
        symbols:
        - is
        - bwtsw
    inputBinding:
      prefix: "-a"

arguments:
  - $("-p" + inputs.InputFile.nameroot)

baseCommand: [bwa, index]

outputs:

  index:
    type: File
    secondaryFiles: 
      - ^.amb
      - ^.ann
      - ^.pac
      - ^.sa
    outputBinding:
      glob: $(inputs.InputFile.nameroot + '.bwt')

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
