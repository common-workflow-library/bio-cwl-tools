#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "sinaiiidgst/lancet:latest"
  InlineJavascriptRequirement: {}

inputs:
  # REQUIRED ARGS

  TumorInput:
    type: File
    inputBinding:
      prefix: "--tumor"
    secondaryFiles:
      - .bai

  NormalInput:
    type: File
    inputBinding:
      prefix: "--normal"
    secondaryFiles:
      - .bai

  Reference: 
    type: File
    format: "http://edamontology.org/format_1929"
    inputBinding:
      prefix: "--ref"

  GenomicRegion:
    type:
      - "null"
      - type: record
        name: GenomicRegion
        fields:
          Chromosome:
            type: string
          RegionStart:
            type: int
          RegionEnd:
            type: int
          GenomicRegion:
            type: string
         # default: $(inputs.Chromosome + ":" + inputs.RegionStart + "-" + inputs.RegionEnd)
            inputBinding:
              prefix: "--reg"
              valueFrom: $(inputs.GenomicRegion.Chromosome + ":" + inputs.GenomicRegion.RegionStart + "-" + inputs.GenomicRegion.RegionEnd)
    
  BedFile:
    type: File?
    inputBinding:
      prefix: "--bed"

# OPTIONAL ARGS

  MinKmer:
    type: int?
    inputBinding:
      prefix: "--min-k"

  MaxKmer:
    type: int?
    inputBinding:
      prefix: "--max-k"

  TrimLowQuality:
    type: int?
    inputBinding:
      prefix: "--trim-lowqual"

  MinBaseQual:
    type: int?
    inputBinding:
      prefix: "--min-base-qual"

  QualityRange:
    type: string?
    inputBinding:
      prefix: "--quality-range"

  MinMapQual:
    type: int?
    inputBinding:
      prefix: "--min-map-qual"

  ASXSDifMax:
    type: int?
    inputBinding:
      prefix: "--max-as-xs-diff"

  MaxTipLen:
    type: int?
    inputBinding:
      prefix: "--tip-len"

  MinCovThres:
    type: int?
    inputBinding:
      prefix: "--cov-thr"

  MinCovRatio:
    type: float?
    inputBinding:
      prefix: "--cov-ratio"

  LowCovThres:
    type: int?
    inputBinding:
      prefix: "--low-cov"

  MaxAvgCov:
    type: int?
    inputBinding:
      prefix: "--max-avg-cov"

  WinSize:
    type: int?
    inputBinding:
      prefix: "--window-size"

  Padding:
    type: int?
    inputBinding:
      prefix: "--padding"

  DFSLimit:
    type: int?
    inputBinding:
      prefix: "--dfs-limit"

  MaxIndelLen:
    type: int?
    inputBinding:
      prefix: "--max-inde-len"

  MaxMismatch:
    type: int?
    inputBinding:
      prefix: "--max-mismatch"

  Threads:
    type: int?
    inputBinding:
      prefix: "--num-threads"

  NodeStrLen:
    type: int?
    inputBinding:
      prefix: "--node-str-len"

  MinAltCountTumor:
    type: int?
    inputBinding:
      prefix: "--min-alt-count-tumor"

  MaxAltCountNormal:
    type: int?
    inputBinding:
      prefix: "--max-alt-count-normal"

  MinVAFTumor:
    type: float?
    inputBinding:
      prefix: "--min-vaf-tumor"

  MaxVAFNormal:
    type: float?
    inputBinding:
      prefix: "--max-vaf-normal"

  MinCovTumor:
    type: int?
    inputBinding:
      prefix: "--min-coverage-tumor"

  MaxCovTumor:
    type: int?
    inputBinding:
      prefix: "--max-coverage-tumor"

  MinCovNormal:
    type: int?
    inputBinding:
      prefix: "--min-coverage-normal"

  MaxCovNormal:
    type: int?
    inputBinding:
      prefix: "--max-coverage-normal"

  MinPhredFisher:
    type: float?
    inputBinding:
      prefix: "--min-phred-fisher"

  MinPhredFisherSTR:
    type: float?
    inputBinding:
      prefix: "--min-phred-fisher-str"

  MinStrandBias:
    type: float?
    inputBinding:
      prefix: "--min-strand-bias"

  MaxUnitLen:
    type: int?
    inputBinding:
      prefix: "--max-unit-length"

  MinReportUnit:
    type: int?
    inputBinding:
      prefix: "--min-report-unit"

  MinReportLen:
    type: int?
    inputBinding:
      prefix: "--min-report-len"

  DistFrSTR:
    type: int?
    inputBinding:
      prefix: "--dist-from-str"

  ActRegOff:
    type: boolean?
    inputBinding:
      prefix: "--active-region-off"

  KmerRecov:
    type: boolean?
    inputBinding:
      prefix: "--kmer-recovery"

  PrintGraph:
    type: boolean?
    inputBinding:
      prefix: "--print-graph"

  Verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"

baseCommand: ["lancet"]

outputs:
  vcf:
    type: stdout

stdout: "lancet-out.vcf"
