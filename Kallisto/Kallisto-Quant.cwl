#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0"

inputs:
  InputReads:
    type: File[]
    format: http://edamontology.org/format_1930 # FASTA
    inputBinding:
      position: 200

  Index:
    type: File
    inputBinding:
      position: 1
      prefix: "--index"

  isSingle:
    type: boolean
    inputBinding:
      position: 2
      prefix: "--single"

  #Optional Inputs

  isBias:
    type: boolean?
    inputBinding:
      prefix: "--bias"

  isFusion:
    type: boolean?
    inputBinding:
      prefix: "--fusion"

  isSingleOverhang:
    type: boolean?
    inputBinding:
      prefix: "--single-overhang"
  
  FragmentLength:
    type: double?
    inputBinding:
      separate: false
      prefix: "--fragment-length="
  
  StandardDeviation:
    type: double?
    inputBinding:
      prefix: "--sd"
  
  BootstrapSamples:
    type: int?
    inputBinding:
      separate: false
      prefix: "--bootstrap-samples="
  
  Seed:
    type: int?
    inputBinding:
      prefix: "--seed"

#Using record inputs to create mutually exclusive inputs
  Strand:
    type:
      - "null"
      - type: record
        name: forward
        fields:
          forward:
              type: boolean
              inputBinding:
                prefix: "--fr-stranded"

      - type: record
        name: reverse
        fields:
          reverse:
            type: boolean
            inputBinding:
              prefix: "--rf-stranded"

  PseudoBam:
    type: boolean?
    inputBinding:
      prefix: "--pseudobam"

#Using record inputs to create dependent inputs
  
  GenomeBam:
    type:
      - "null"
      - type: record
        name: genome_bam
        fields:
          genomebam:
            type: boolean
            inputBinding:
              prefix: "--genomebam"

          gtf:
            type: File
            inputBinding:
              prefix: "--gtf"

          chromosomes:
            type: File
            inputBinding:
              prefix: "--chromosomes"

baseCommand: [ kallisto, quant ]

arguments: [ "--output-dir", out ]

outputs:

  quantification_h5:
    type: File
    outputBinding:
      glob: out/abundances.h5

# Long form method for defining optional outputs

  quantification_tsv:
    type: File
    outputBinding:
      glob: out/abundances.tsv

  bam:
    type: ["null", File]
    outputBinding:
      glob: "out/*.bam"

  fusions:
    type: ["null", File]
    outputBinding:
      glob: "fusion.txt"
  
