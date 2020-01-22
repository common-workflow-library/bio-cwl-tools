#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/pizzly:0.37.3--0"

inputs:
  InputFile:
    type: File
    inputBinding:
      position: 100

  Reference:
    type: File
    inputBinding:
      prefix: "--fasta"

  Kmer:
    type: int
    inputBinding:
      prefix: "-k"

  Gtf:
    type: File
    inputBinding:
      prefix: "--gtf"

  Output:
    type: string
    default: "pizzly_out"
    inputBinding:
      prefix: "--output"
      valueFrom: "pizzly_out"

  #Optional Inputs

  Cache:
    type: string?
    default: "index.cache.txt"
    inputBinding:
      prefix: "--cache"
  
  InsertSize:
    type: int?
    inputBinding:
      prefix: "--insert-size"

  isIgnoreProtein:
    type: boolean?
    inputBinding:
      prefix: "--ignore-protein"

baseCommand: ["pizzly"]

outputs:

  Fusion_fasta:
    type: File
    outputBinding:
      glob: .fusions.fasta

  Fusion_json:
    type: File
    outputBinding:
      glob: .json
