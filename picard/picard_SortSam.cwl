#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/picard:2.22.2--0
requirements:
  InlineJavascriptRequirement: {}

baseCommand: [ picard, SortSam ]

arguments:
  - prefix: OUTPUT=
    separate: false
    valueFrom: |
      ${ if(inputs.sort_order == "coordinate") { return (inputs.inputFile.nameroot)+".bam";} else { return (inputs.inputFile.nameroot)+".sam"; } }

inputs:
  alignments:
    type: File
    inputBinding:
      prefix: INPUT=
      separate: false

  sort_order:
    type:
      - 'null'
      - type: enum
        symbols:
          - queryname
          - coordinate
          - duplicate
    default: coordinate
    doc: 'coordinate (bam) or queryname (sam)'
    inputBinding:
      prefix: SORT_ORDER=
      separate: false

  validation_stringency:
    default: LENIENT
    doc: Validation stringency for all SAM files read by this program.  Setting stringency
      to SILENT can improve performance when processing a BAM file in which variable-length
      data (read, qualities, tags) do not otherwise need to be decoded.
    type:
    - 'null'
    - type: enum
      symbols:
      - STRICT
      - LENIENT
      - SILENT
    inputBinding:
      prefix: VALIDATION_STRINGENCY=
      separate: false

outputs:
  sorted_alignments:
    type: File
    outputBinding:
      glob: '*.*am'
