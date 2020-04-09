class: CommandLineTool
cwlVersion: v1.0
id: picard_MarkDuplicates
doc: "examine aligned records in BAM datasets to locate duplicate molecules"

requirements:
  InlineJavascriptRequirement: {}
hints:
  - class: DockerRequirement
    dockerPull: broadinstitute/picard:latest

baseCommand: [java, -jar, /usr/picard/picard.jar]
arguments:
  - MarkDuplicates
#  - VERBOSY=ERROR

inputs:
  - id: inputFile
    type: File
    inputBinding:
      position: 1
      prefix: INPUT=
      separate: false

  - id: outFile_name
    type: string?
    inputBinding:
      position: 2
      prefix: OUTPUT=
      valueFrom: ${return (inputs.inputFile.nameroot)+"_markduplicates.bam"}
      separate: false

  - id: metric_file_name
    type: string?
    default: metrics_file.txt
    inputBinding:
      position: 3
      prefix: METRICS_FILE=
      separate: false

  - id: comment
    type: string[]?
    inputBinding:
      position: 4
      prefix: COMMENT=
      separate: false

  - id: remove_duplicates
    type: string?
    default: 'false'
    inputBinding:
      position: 5
      prefix: REMOVE_DUPLICATES=
      separate: false

  - id: assume_sorted
    type: string?
    default: 'true'
    inputBinding:
      position: 6
      prefix: ASSUME_SORTED=
      separate: false

  - id: duplicate_scoring_strategy
    type: string?
    default: 'SUM_OF_BASE_QUALITIES'
    doc: 'SUM_OF_BASE_QUALITIES OR TOTAL_MAPPED_REFERENCE_LENGH'
    inputBinding:
      position: 7
      prefix: DUPLICATE_SCORING_STRATEGY=
      separate: false

  - id: reg_name_regex
    type: string?
    inputBinding:
      position: 8
      prefix: READ_NAME_REGEX=
      separate: false

  - id: optical_duplicate_pixel_distance
    type: int?
    default: 100
    doc: '(0;500)'
    inputBinding:
      position: 9
      prefix: OPTICAL_DUPLICATE_PIXEL_DISTANCE=
      separate: false

  - id: barcode_tag
    type: string?
    inputBinding:
      position: 10
      prefix: BARCODE_TAG=
      separate: false

  - id: validation_stringency
    type: string?
    default: LENIENT
    doc: 'LENIENT, SILENT or STRICT'
    inputBinding:
      position: 11
      prefix: VALIDATION_STRINGENCY=
      separate: false


outputs:
  - id: outFile
    type: File
    outputBinding:
      glob: '*.bam'
  - id: metrics_file
    type: File
    outputBinding:
      glob: '*.txt'
