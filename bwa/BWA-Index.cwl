#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa:0.7.17--ha92aebf_3
  SoftwareRequirement:
    packages:
      bwa:
        version: [ "0.7.17" ]
        specs: [ https://identifiers.org/biotools/bwa ]
          
inputs:
  sequences:
    type: File
    format: edam:format_1929  # FASTA
    inputBinding:
      position: 200

#Optional arguments

  algo_type:
    type: 
      - "null"
      - type: enum
        symbols:
        - is
        - bwtsw
    inputBinding:
      prefix: "-a"
  index_name:
    type: string?

# the expressions for the index name prefix are needed because 'default:' on an optional
# parameter does not accept an expression
arguments:
  - -p
  - '$(((inputs.index_name !== null) ? inputs.index_name : inputs.sequences.nameroot))'

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
      glob: "*.bwt"

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
