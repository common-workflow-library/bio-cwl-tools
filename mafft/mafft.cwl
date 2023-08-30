#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: Mafft
doc: |-
  MAFFT (Multiple Alignment using Fast Fourier Transform) is a high speed multiple sequence alignment program.
$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/

inputs:
  add:
    doc: add unaligned full-length sequences into an existing alignment
    type: File?
    format: edam:format_1929
    inputBinding:
      prefix: --add
  addfragments:
    doc: add unaligned fragmentary sequences into an existing alignment
    type: File?
    format: edam:format_1929
    inputBinding:
      prefix: --addfragments
  anysymbol:
    doc: allow symbols not part of the standard IUPAC nucleotide or protein alphabets
    type: boolean?
    inputBinding:
      prefix: --anysymbol
  auto:
    doc: auto-select alignment strategy
    type: boolean?
    default: true
    inputBinding:
      prefix: --auto
  no_save_memory:
    doc: always apply normal DP even for long alingments
    type: boolean?
    inputBinding:
      prefix: --nomemsave
  save_memory:
    doc: use linear-space DP algorithm (approximately two times slower than normal
      DP)
    type: boolean?
    inputBinding:
      prefix: --memsave
  sequences:
    label: Sequences to align
    type: File
    format: edam:format_1929
    inputBinding:
      position: 1

outputs:
  alignment:
    type: File
    format: edam:format_1929
    outputBinding:
      glob: $(inputs.sequences.nameroot).alignment.fasta
    streamable: true
stdout: $(inputs.sequences.nameroot).alignment.fasta

baseCommand: mafft
arguments:
- prefix: --thread
  valueFrom: $(runtime.cores)

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/mafft:7.458--h516909a_0
  ResourceRequirement:
    coresMin: 8
    ramMin: 40000
  SoftwareRequirement:
    packages:
      mafft:
        specs:
          - https://identifiers.org/biotools/MAFFT
          - https://anaconda.org/bioconda/mafft
        version: [ "7.458" ]
$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf
- https://edamontology.org/EDAM_1.18.owl
s:author:
- class: s:Person
  s:identifier: https://orcid.org/0000-0002-2961-9670
  s:name: Michael R. Crusoe
- class: s:Person
  s:email: mailto:pvh@sanbi.ac.za
  s:identifier: https://orcid.org/0000-0001-6553-5274
  s:name: Peter van Heusden
