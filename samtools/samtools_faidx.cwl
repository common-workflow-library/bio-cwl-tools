#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.2-0
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]

baseCommand: [ samtools, faidx ]

inputs:
  sequences:
    type: File
    doc: Input FASTA file
    format: edam:format_1929

arguments:
   - $(inputs.sequences.basename)

outputs:
  sequences_with_index:
    type: File
    format: $(inputs.sequences.format)
    secondaryFiles:
     - .fai
    outputBinding:
      glob: $(inputs.sequences.basename)
  sequences_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename).fai

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
