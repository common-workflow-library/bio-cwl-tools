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
      prefix: "--input"

  Output: 
    type: string
    default: "ApplyBQSR.bam"
    inputBinding:
      prefix: "--output" 
      valueFrom: "ApplyBQSR.bam"

  BaseRecalFile: 
    type: File
    inputBinding:
      prefix: "--bqsr-recal-file" 

  # OPTIONAL ARGS

  isOutputSamProgramRecord:
    type: boolean?
    inputBinding:
      prefix: "--add-output-sam-program-record"

  isOutputVCFCommandLine:
    type: boolean?
    inputBinding:
      prefix: "--add-output-vcf-command-line"

  ArgumentsFile:
    type: File?
    inputBinding:
      prefix: "--arguments_file"

  CloudIndexPrefetchBuilder:
    type: int?
    inputBinding:
      prefix: "--cloud-index-prefetch-buffer"

  CloudPrefetchBuffer:
    type: int?
    inputBinding:
      prefix: "--cloud-prefetch-buffer"

  isCreateBamIndex:
    type: boolean?
    inputBinding:
      prefix: "--create-output-bam-index"

  isCreateMD5:
    type: boolean?
    inputBinding:
      prefix: "--create-output-bam-md5"

  isCreateOutputVariantIndex:
    type: boolean?
    inputBinding:
      prefix: "--create-output-variant-index"

  isCreateOutputVariantMD5:
    type: boolean?
    inputBinding:
      prefix: "--create-output-variant-md5"

  isDisableBamIndexCaching:
    type: boolean?
    inputBinding:
      prefix: "--disable-bam-index-caching"

  isDisableReadFilter:
    type: boolean?
    inputBinding:
      prefix: "--disable-read-filter"

  isDisableSequenceDictionaryValidation:
    type: boolean?
    inputBinding:
      prefix: "--disable-sequence-dictionary-validation"

  isDisableToolDefaultReadFilters:
    type: boolean?
    inputBinding:
      prefix: "--disable-tool-default-read-filters"

  isEmitOriginalQuals:
    type: boolean?
    inputBinding:
      prefix: "--emit-original-quals"

  ExcludeIntervals:
    type: 
      - "null"
      - type: array
        items: string
        inputBinding:
          prefix: "--exclude-intervals"
      - File
    inputBinding:
      prefix: "--exclude-intervals"

  GatkConfigFile:
    type: File?
    inputBinding:
      prefix: "--gatk-config-file"

  GCSMaxRetries:
    type: int?
    inputBinding:
      prefix: "--gcs-max-retries"
  
  GlobalScopePrior:
    type: double?
    inputBinding:
      prefix: "--global-qscore-prior"

  IntervalExclusionPadding:
    type: int?
    inputBinding:
      prefix: "--interval-exclusion-padding"

  IntervalMergingRule:
    type: 
      - "null"
      - type: enum
        symbols:
          - ALL
          - OVERLAPPING_ONLY
    inputBinding:
      prefix: "--interval-merging-rule"

  IntervalPadding:
    type: int?
    inputBinding:
      prefix: "--interval-padding"

  IntervalSetRule:
    type: 
      - "null"
      - type: enum
        symbols:
          - UNION
          - INTERSECTION

  Intervals:
    type: string[]?
    inputBinding:
      prefix: "--intervals"

  JavaOptions: 
    type: string?
    inputBinding:
      prefix: "--java_options" 
      position: -2

  isLenient:
    type: boolean?
    inputBinding:
      prefix: "--lenient"

  PreserveQscoresLessThan:
    type: int?
    inputBinding:
      prefix: "--preserve-qscores-less-than"

  QuantizeQuals:
    type:
      - "null"
      - type: record
        name: QuantizeQuals
        fields:
          QuantizeQuals:
              type: int
              inputBinding:
                prefix: "--quantize-quals"
      - type: record
        name: RoundDownQuantize
        fields:
          RoundDownQuantize:
            type: boolean?
            inputBinding:
              prefix: "--round-down-quantized"

          StaticQuantizedQuals:
            type: boolean?
            inputBinding:
              prefix: "--static-quantized-quals"

  ReadFilter:
    type: string[]?
    inputBinding:
      prefix: "--read-filter"
      itemSeparator: ","

  ReadIndex:
    type: string[]?
    inputBinding:
      prefix: "--read-index"
      itemSeparator: ","

  ReadValidationStringency:
    type: 
      - "null"
      - type: enum
        symbols:
          - STRICT
          - LENIENT
          - SILENT
    inputBinding:
      prefix: "--read-validation-stringency"

  Reference:
    type: string?
    inputBinding:
      prefix: "--reference"

  SecondsBetweenProgUpdates:
    type: int?
    inputBinding:
      prefix: "--seconds-between-progress-updates"

  SequenceDictionary:
    type: File?
    inputBinding:
      prefix: "--sequence-dictionary"

  isUseOriginalQualities:
    type:  boolean?
    inputBinding:
      prefix: "--use-original-qualities"

  isJDKDeflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_DEFLATER"  
  
  isJDKInflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_INFLATER"

baseCommand: ["/gatk/gatk"]

arguments: 
  - valueFrom: "ApplyBQSR"
    position: -1

outputs:
  alignment:
    type: File
    outputBinding:
      glob: "ApplyBQSR.bam"

  index:
    type: ["null", File]
    outputBinding:
      glob: "ApplyBQSR.bai"
  vcf:
    type: ["null", File]
    outputBinding:
      glob: "*.vcf"

  

