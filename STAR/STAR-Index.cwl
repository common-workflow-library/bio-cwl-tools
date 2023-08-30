#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/star:2.7.5c--0
  SoftwareRequirement:
    packages:
      star:
        specs: [ https://identifiers.org/biotools/star ]
        version: [ "2.7.5c" ]

inputs:

  InputFiles:
    format: edam:format_1929  # FASTA
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

$namespaces:
  edam: http://edamontology.org/
$schemas: [ https://edamontology.org/EDAM_1.25.owl ]
