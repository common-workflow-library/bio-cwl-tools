#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: docker pull quay.io/biocontainers/samtools:1.2-0
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]

baseCommand: [ samtools, faidx ]

inputs:
  sequences:
    type: File
    doc: Input FASTA file

arguments:
   - $(inputs.sequences.basename)

outputs:
  sequences_with_index:
    type: File
    secondaryFiles:
     - .fai
    outputBinding:
      glob: $(inputs.sequences)
  sequences_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename).fai
