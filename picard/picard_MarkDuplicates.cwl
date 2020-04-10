#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Removal of duplicates from aligned reads.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/picard:2.22.2--0
  
baseCommand: [ picard, MarkDuplicates ]

arguments:
  - OUTPUT=$(inputs.alignments.nameroot)_markduplicates$(inputs.alignments.nameext)
  - METRICS_FILE=$(inputs.alignments.nameroot)_markduplicates.metrics

stderr: $(inputs.alignments.nameroot).markduplicates.log

inputs:
  alignments:
    doc: SAM or BAM format alignment file
    format:
      - edam:format_2573  # SAM
      - edam:format_2572  # BAM
    type: File
    inputBinding:
      prefix: "INPUT="
      separate: false

  alignments_are_sorted:
    type: boolean
    inputBinding:
      prefix: ASSUME_SORTED=TRUE

  remove_duplicates:
    doc: |
     If true do not write duplicates to the output file instead of writing them
     with appropriate flags set.
    type: boolean
    inputBinding:
      prefix: REMOVE_DUPLICATES=TRUE

  validation_stringency:
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

  comment:
    doc: Comment(s) to include in the output file's header
    type:
     - 'null'
     - type: array
       items: string
       inputBinding:
         prefix: COMMENT=
         separate: false

  duplicate_scoring_strategy:
    type:
      - 'null'
      - type: enum
        symbols:
          - SUM_OF_BASE_QUALITIES
          - TOTAL_MAPPED_REFERENCE_LENGTH
          - RANDOM
    inputBinding:
      prefix: DUPLICATE_SCORING_STRATEGY=
      separate: false

  read_name_regex:
    type: string?
    inputBinding:
      prefix: READ_NAME_REGEX=
      separate: false

  optical_duplicate_pixel_distance:
    type: int?
    default: 100
    doc: '(0;500)'
    inputBinding:
      prefix: OPTICAL_DUPLICATE_PIXEL_DISTANCE=
      separate: false

  barcode_tag:
    type: string?
    inputBinding:
      prefix: BARCODE_TAG=
      separate: false

outputs:
  alignments:
    type: File
    format: $(inputs.alignments.format)
    outputBinding:
      glob: $(inputs.alignments.nameroot)_markduplicates$(inputs.alignments.nameext)
  log:
    type: stderr
  metrics:
    type: File
    outputBinding:
      glob: $(inputs.alignments.nameroot)_markduplicates.metrics

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
