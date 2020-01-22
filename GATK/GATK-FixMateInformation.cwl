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
    default: "fixmate.dat"
    inputBinding:
      prefix: "--OUTPUT" 
      valueFrom: $("fixmate.dat")

  # OPTIONAL ARGS

  isAddMateCigar:
    type: boolean?
    inputBinding:
      prefix: "--ADD_MATE_CIGAR"

  ArgumentsFile:
    type: File?
    inputBinding:
      prefix: "--arguments_file"

  isAssumedSorted:
    type: boolean?
    inputBinding:
      prefix: "--ASSUME_SORTED"

  isIgnoreMissingMates:
    type: boolean?
    inputBinding:
      prefix: "--IGNORE_MISSING_MATES"

  SortOrder:
    type: 
      - "null"
      - type: enum
        symbols:
          - unsorted
          - queryname
          - coordinate
          - duplicate
          - unknown
    inputBinding:
      prefix: "--SORT_ORDER"

  JavaOptions: 
    type: string?
    inputBinding:
      prefix: "--java_options" 
      position: -2

  isJDKDeflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_DEFLATER"  
  
  isJDKInflator:
    type: boolean?
    inputBinding:
      prefix: "--USE_JDK_INFLATER"

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

baseCommand: ["/gatk/gatk"]

arguments: 
  - valueFrom: "FixMateInformation"
    position: -1

outputs:
  alignment:
    type: File
    outputBinding:
      glob: fixmate.dat
  index:
    type: ["null", File]
    outputBinding:
      glob: ^.bai
  MD5:
    type: ["null", File]
    outputBinding:
      glob: ^.MD5