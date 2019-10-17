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
      prefix: "-I"

  Reference:
    type: File
    inputBinding:
      prefix: "-R"
    secondaryFiles:
      - .dict

  Output:
    type: string
    default: "recal_data.table"
    inputBinding:
      prefix: "--output"
      valueFrom: "recal_data.table"

  KnownSites: 
    type: File
    inputBinding:
      prefix: "--known-sites" 
    secondaryFiles:
      - .idx

  # OPTIONAL ARGS

  Covariates:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
        prefix: "--covariate"

  IndelsContextSize:
    type: int?
    inputBinding:
      prefix: "--indels_context_size"

  MaximumCycleValue:
    type: int?
    inputBinding:
      prefix: "--maximum_cycle_value"

  MismatchesContextSize:
    type: int?
    inputBinding:
      prefix: "--mismatches_context_size"
 
  SolidNoCallStrategy:
    type: 
      - "null"
      - type: enum
        symbols:
          - THROW_EXCEPTION
          - LEAVE_READ_UNRECALIBRATED
          - PURGE_READ
    inputBinding:
      prefix: "--solid_nocall_strategy"

  SolidRecalMode:
    type:
      - "null"
      - type: enum
        symbols:
          - DO_NOTHING
          - SET_Q_ZERO
          - SET_Q_ZERO_BASE_N
          - REMOVE_REF_BIAS
    inputBinding:
      prefix: "--solid_recal_mode"

  isListCovariates:
    type: boolean?
    inputBinding:
      prefix: "--list"

  isLowMemoryMode:
    type: boolean?
    inputBinding:
      prefix: "--lowMemoryMode"

  #No standard covariates could be an issue

  isSortAllColumns:
    type: boolean?
    inputBinding:
      prefix: "--sort_by_all_columns"

  BinaryTagName:
    type: string?
    inputBinding:
      prefix: "--binary_tag_name"

  bqsrBAQGapOpenPenalty:
    type: float?
    inputBinding:
      prefix: "--bqsrBAQGapOpenPenalty"

  DeletionsDefaultQuality:
    type: int?
    inputBinding:
      prefix: "--deletions_default_quality"

  InsertionsDefaultQuality:
    type: int?
    inputBinding:
      prefix: "--insertions_default_quality"

  LowQualityTail:
    type: int?
    inputBinding:
      prefix: "--low_quality_tail"

  MismatchesDefaultQuality:
    type: int?
    inputBinding:
      prefix: "--mismatches_default_quality"

  QuantizingLevels:
    type: int?
    inputBinding:
      prefix: "--quantizing_levels"

  RunWithoutdbSNP:
    type: boolean?
    inputBinding:
      prefix: "--run_without_dbsnp_potentially_ruining_quality"
  

arguments: 
  - valueFrom: "BaseRecalibrator"
    position: -1

baseCommand: ["/gatk/gatk"]

outputs:
  table:
    type: File
    outputBinding:
      glob: $(inputs.Output)
    