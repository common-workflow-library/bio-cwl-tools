#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "broadinstitute/gatk:4.1.1.0"
  InlineJavascriptRequirement: {}

inputs:
  # REQUIRED ARGS

  InputFile:
    type: File
    inputBinding:
      prefix: "--INPUT"

  Output: 
    type: string
    default: $("MarkDuplicatesOut" + inputs.InputFile.nameext)
    inputBinding:
      prefix: "--OUTPUT" 
      valueFrom: $("MarkDuplicatesOut" + inputs.InputFile.nameext)

  MetricsFile: 
    type: string
    default: "MarkDuplicatesMetrics.txt"
    inputBinding:
      prefix: "--METRICS_FILE" 
      valueFrom: "MarkDuplicatesMetrics.txt"

  # OPTIONAL ARGS

  ArgumentsFile:
    type: File[]?
    inputBinding:
      prefix: "--arguments_file"

  isAssumeSortOrder:
    type: boolean?
    inputBinding:
      prefix: "--ASSUME_SORT_ORDER"

  isBarcodeTag:
    type: boolean?
    inputBinding:
      prefix: "--BARCODE_TAG"
  
  isClearDT:
    type: boolean?
    inputBinding:
      prefix: "--CLEAR_DT"

  Comment:
    type: string[]?
    inputBinding:
      itemSeparator: ","
      prefix: "--COMMENT"
  
  DuplicateScoringStrategy:
    type: 
      - "null"
      - type: enum
        symbols:
          - SUM_OF_BASE_QUALITIES
          - TOTAL_MAPPED_REFERENCE_LENGTH
          - RANDOM
    inputBinding:
      prefix: "--DUPLICATE_SCORING_STRATEGY"

  JavaOptions: 
    type: string?
    inputBinding:
      prefix: "--java_options" 
      position: -2

  MaxFileHandles:
    type: int?
    inputBinding:
      prefix: "--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP"

  MaxOpticalDuplicate:
    type: int?
    inputBinding:
      prefix: "--MAX_OPTICAL_DUPLICATE_SET_SIZE"

  MaxSequenceDiskReadMap:
    type: int?
    inputBinding:
      prefix: "--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP"

  OpticalDuplicatePixelDistance:
    type: int?
    inputBinding:
      prefix: "--OPTICAL_DUPLICATE_PIXEL_DISTANCE"

  ProgramGroupCommandLine:
    type: string?
    inputBinding:
      prefix: "--PROGRAM_GROUP_COMMAND_LINE"

  ProgramGroupName:
    type: string?
    inputBinding:
      prefix: "--PROGRAM_GROUP_NAME"

  ProgramGroupVersion:
    type: string?
    inputBinding:
      prefix: "--PROGRAM_GROUP_VERSION"

  ProgramRecordId:
    type: string?
    inputBinding:
      prefix: "--PROGRAM_RECORD_ID"

  ReadNameRegex:
    type: string?
    inputBinding:
      prefix: "--READ_NAME_REGEX"

  ReadOneBarcodeTag:
    type: string?
    inputBinding:
      prefix: "--READ_ONE_BARCODE_TAG"

  ReadTwoBarcodeTag:
    type: string?
    inputBinding:
      prefix: "--READ_TWO_BARCODE_TAG"

  isRemoveDuplicates:
    type: boolean?
    inputBinding:
      prefix: "--REMOVE_DUPLICATES"

  isRemoveSequenceDuplicates:
    type: boolean?
    inputBinding:
      prefix: "--REMOVE_SEQUENCING_DUPLICATES"

  SortingCollectionSizeRatio:
    type: boolean?
    inputBinding:
      prefix: "--SORTING_COLLECTION_SIZE_RATIO"

  isTagDuplicateSetMembers:
    type: boolean?
    inputBinding:
      prefix: "--TAG_DUPLICATE_SET_MEMBERS"

  TaggingPolicy:
    type: 
      - "null"
      - type: enum
        symbols:
          - DontTag
          - OpticalOnly
          - All
    inputBinding:
      prefix: "--TAGGING_POLICY"

  AddPGTagToReads:
    type: boolean?
    inputBinding:
      prefix: "--ADD_PG_TAG_TO_READS"

  CompresionLevel:
    type: int?
    inputBinding:
      prefix: "--COMPRESSION_LEVEL"

  isCreateIndex:
    type: boolean?
    inputBinding:
      prefix: "--CREATE_INDEX"

  isCreateMD5:
    type: boolean?
    inputBinding:
      prefix: "--CREATE_MD5_FILE"

  GH4Secrets:
    type: File?
    inputBinding:
      prefix: "--GA4GH_CLIENT_SECRETS"

  MaxRecordsRam:
    type: int?
    inputBinding:
      prefix: "--MAX_RECORDS_IN_RAM"

  isQuiet:
    type: boolean?
    inputBinding:
      prefix: "--QUIET"

  ReferenceSequence:
    type: File?
    inputBinding:
      prefix: "--REFERENCE_SEQUENCE"

  isJDKDeflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_DEFLATER"

  isJDKInflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_INFLATER"

  ValidationStringency: 
    type: 
     - "null"
     - type: enum
       symbols:
        - STRICT
        - LENIENT
        - SILENT
    inputBinding:
      prefix: "--VALIDATION_STRINGENCY"

  Verbosity: 
    type: 
     - "null"
     - type: enum
       symbols:
        - ERROR
        - WARNING
        - INFO
        - DEBUG
    inputBinding:
      prefix: "--VERBOSITY"
 
  Version: 
    type: boolean?
    inputBinding:
      prefix: "--version" 

arguments: 
  - valueFrom: "MarkDuplicates"
    position: -1

baseCommand: ["/gatk/gatk"]

outputs:
  alignment:
    type: File
    outputBinding:
      glob: $("MarkDuplicatesOut" + inputs.InputFile.nameext)
  metrics:
    type: File
    outputBinding:
      glob: "MarkDuplicatesMetrics.txt"
  index:
    type: ["null", File]
    outputBinding:
      glob: "*.bai"
  vcf:
    type: ["null", File]
    outputBinding:
      glob: "*.vcf"

