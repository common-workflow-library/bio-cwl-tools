#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/star:2.7.5c--0"

inputs:

  InputFiles:
    format: http://edamontology.org/format_1930
    type: File[]
    inputBinding:
      prefix: "--genomeFastaFiles"

  IndexName:
    type: string
    inputBinding:
      prefix: "--genomeDir"
      valueFrom: ./$(self)

#Optional Inputs

  Gtf:
    type: File?
    inputBinding:
      prefix: "--sjdbGTFfile"

  Overhang:
    type: int?
    inputBinding:
      prefix: "--sjdbOverhang"

  Junctions:
    type: File?
    inputBinding:
      prefix: "--sjdbFileChrStartEnd"

  GenomeSize:
    type: int?
    inputBinding:
      prefix: "--genomeSAindexNbases"

  GenomeBits:
    type: int?
    inputBinding:
      prefix: "--genomeChrBinNbits"

baseCommand: [STAR, --runMode, genomeGenerate]

arguments: [--runThreadN, $(runtime.cores)]

outputs:
  indexes:
    type: Directory
    outputBinding:
      glob: ./$(inputs.IndexName)/
