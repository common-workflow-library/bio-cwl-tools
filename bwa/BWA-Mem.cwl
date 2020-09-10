#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
  StepInputExpressionRequirement: {}

inputs:
  InputFile:
    type: File[]
    format:
      - edam:format_1930 # FASTQ (no quality score encoding specified)
      - edam:format_1931 # FASTQ-Illumina
      - edam:format_1932 # FASTQ-Sanger
      - edam:format_1933 # FASTQ-Solexa
    inputBinding:
      position: 201
    
  Index:
    type: File
    secondaryFiles:
      - ^.amb
      - ^.ann
      - ^.pac
      - ^.sa
    inputBinding:
      position: 200
      valueFrom: $(self.dirname + '/' + self.nameroot)
    
#Optional arguments

  Threads:
    type: int?
    inputBinding:
      prefix: "-t"

  MinSeedLen:
    type: int?
    inputBinding:
      prefix: "-k"
  
  BandWidth:
    type: int?
    inputBinding:
      prefix: "-w"

  ZDropoff:
    type: int?
    inputBinding:
      prefix: "-d"

  SeedSplitRatio:
    type: float?
    inputBinding:
      prefix: "-r"
    
  MaxOcc:
    type: int?
    inputBinding:
      prefix: "-c"

  MatchScore:
    type: int?
    inputBinding:
      prefix: "-A"

  MmPenalty:
    type: int?
    inputBinding:
      prefix: "-B"

  GapOpenPen:
    type: int?
    inputBinding:
      prefix: "-O"

  GapExtPen:
    type: int?
    inputBinding:
      prefix: "-E"

  ClipPen:
    type: int?
    inputBinding:
      prefix: "-L"

  UnpairPen:
    type: int?
    inputBinding:
      prefix: "-U"

  RgLine:
    type: string?
    inputBinding:
      prefix: "-R"

  VerboseLevel:
    type: int?
    inputBinding:
      prefix: "-v"

  isOutSecAlign:
    type: boolean?
    inputBinding:
      prefix: "-a"

  isMarkShortSplit:
    type: boolean?
    inputBinding:
      prefix: "-M"

  isUseHardClip:
    type: boolean?
    inputBinding:
      prefix: "-H"

  isMultiplexedPair:
    type: boolean?
    inputBinding:
      prefix: "-p"
      

baseCommand: [bwa, mem]
    
stdout: unsorted_reads.sam

outputs:
  reads_stdout:
    type: stdout
    format: edam:format_2573 # SAM format
    
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl