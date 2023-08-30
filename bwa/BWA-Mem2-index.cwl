#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

inputs:
  sequences:
    type: File
    format: edam:format_1929  # FASTA
 
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_2
  SoftwareRequirement:
    packages:
      bwa-mem2:
        version: [ 2.2.1 ]                                                      
        specs: [ https://identifiers.org/biotools/bwa-mem2 ]

requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.sequences)  
   
baseCommand: [bwa-mem2, index]

arguments:
  - $(inputs.sequences)

outputs:
  indexed_sequences:
    type: File
    format: edam:format_1929  # FASTA
    outputBinding:
      glob: $(inputs.sequences.basename)
    secondaryFiles:
      - .bwt.2bit.64
      - .ann
      - .amb
      - .pac
      - ".0123"

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
