#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/hisat2:2.0.4--py35_1"

inputs:
  # Required Inputs
  RunThreadN:
    type: int
    inputBinding:
      prefix: "--runThreadN"

  GenomeDir:
    type: Directory
    inputBinding:
      prefix: "--genomeDir"

  ForwardReads:
    format: http://edamontology.org/format_1930
    type:
     - File
     - File[]
    inputBinding:
      prefix: "--readFilesIn"
      separate: false
      itemSeparator: ","
      position: 1
  # If paired-end reads (like Illumina), both 1 and 2 must be provided.
  ReverseReads:
    format: http://edamontology.org/format_1930
    type:
     - "null"
     - File
     - File[]
    inputBinding:
      prefix: ""
      separate: false
      itemSeparator: ","
      position: 2

  # Optional Inputs
  Gtf:
    type: File?
    inputBinding:
      prefix: "--sjdbGTFfile"

  Overhang:
    type: int?
    inputBinding:
      prefix: "--sjdbOverhang"

  OutFilterType:
    type:
     - "null"
     - type: enum
       symbols:
        - Normal
        - BySJout
    inputBinding:
      prefix: "--outFilterType"

  OutFilterIntronMotifs:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - RemoveNoncanonical
        - RemoveNoncanonicalUnannotated
    inputBinding:
      prefix: "--outFilterIntronMotifs"
  
  OutSAMtype:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - "BAM"
        - "BAM Unsorted"
        - "BAM SortedByCoordinate"
        - "BAM Unsorted SortedByCoordinate"
        - "SAM"
        - "SAM Unsorted"
        - "SAM SortedByCoordinate"
        - "SAM Unsorted SortedByCoordinate"
    inputBinding:
      prefix: "--outSAMtype"
  
  ReadFilesCommand:
    type: string?
    inputBinding:
      prefix: "--readFilesCommand"

  AlignIntronMin:
    type: int?
    inputBinding:
      prefix: "--alignIntronMin"
  
  AlignIntronMax:
    type: int?
    inputBinding:
      prefix: "--alignIntronMax"
  
  AlignMatesGapMax:
    type: int?
    inputBinding:
      prefix: "--alignMatesGapMax"

  AlignSJoverhangMin:
    type: int?
    inputBinding:
      prefix: "--alignSJoverhangMin"
  
  AlignSJDBoverhangMin:
    type: int?
    inputBinding:
      prefix: "--alignSJDBoverhangMin"
  
  SeedSearchStartLmax:
    type: int?
    inputBinding:
      prefix: "--seedSearchStartLmax"

  ChimOutType:
    type:
     - "null"
     - type: enum
       symbols:
        - Junctions
        - SeparateSAMold
        - WithinBAM
        - "WithinBAM HardClip"
        - "WithinBAM SoftClip"

  ChimSegmentMin:
    type: int?
    inputBinding:
      prefix: "--chimSegmentMin"
    
  ChimJunctionOverhangMin:
    type: int?
    inputBinding:
      prefix: "--chimJunctionOverhangMin"

  OutFilterMultimapNmax:
    type: int?
    inputBinding:
      prefix: "--outFilterMultimapNmax"
  
  OutFilterMismatchNmax:
    type: int?
    inputBinding:
      prefix: "--outFilterMismatchNmax"

  OutFilterMismatchNoverLmax:
    type: double?
    inputBinding:
      prefix: "--outFilterMismatchNoverLmax"

  OutReadsUnmapped:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - Fastx
    inputBinding:
      prefix: "--outReadsUnmapped"
  
  OutSAMstrandField:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - intronMotif
    inputBinding:
      prefix: "--outSAMstrandField"
  
  OutSAMunmapped:
    type:
     - "null"
     - type: enum
       symbols:
        - None
        - Within
        - "Within KeepPairs"
    inputBinding:
      prefix: "--outSAMunmapped"
  
  OutSAMmapqUnique:
    type: int?
    inputBinding:
      prefix: "--outSAMmapqUnique"
  
  OutSamMode:
    type: 
     - "null"
     - type: enum
       symbols:
        - None
        - Full
        - NoQS
    inputBinding:
      prefix: "--outSAMmode"
  
  LimitOutSAMoneReadBytes:
    type: int?
    inputBinding:
      prefix: "--limitOutSAMoneReadBytes"
  
  OutFileNamePrefix:
    type: string?
    inputBinding:
      prefix: "--outFileNamePrefix"

  GenomeLoad:
    type:
     - "null"
     - type: enum
       symbols:
        - LoadAndKeep
        - LoadAndRemove
        - LoadAndExit
        - Remove
        - NoSharedMemory
    inputBinding:
      prefix: "--genomeLoad"

baseCommand: [STAR, --runmode, alignReads]       

outputs:
  alignment:
    type:
     - File
     - File[]
    outputBinding:
      glob: "*.bam"
  unmapped_reads:
    type: ["null", File]
    outputBinding:
      glob: "Unmapped.out*"
