#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: minimap2

doc: This CWL file defines running minimap2 to align some sequences to a database.
  We assume the database has been indexed. This is not necessary but we will do it
  in our use case

requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/minimap2:2.17--h8b12597_1
  SoftwareRequirement:
    packages:
      minimap2:
        specs:
          - https://identifiers.org/biotools/minimap2
          - https://github.com/lh3/minimap2
        version: [ "2.17" ]
  ResourceRequirement:
    coresMin: 8
    coresMax: 32
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.target.size/(1024*1024*1024) + 20))

inputs:
  preset:
    type:
      - 'null'
      - type: enum
        symbols:
          - map-pb
          - map-ont
          - ava-pb
          - ava-ont
          - asm5
          - asm10
          - asm20
          - splice
          - sr
    inputBinding:
      prefix: "-x"
  outputCIGAR:
    type: boolean?
    label: output CIGAR in PAF
    inputBinding:
      prefix: -c
  miniWinSize:
    type: int?
    label: minimizer window length
    inputBinding:
      prefix: -w
  target:
    type: File
    inputBinding:
      position: 5
  query:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      position: 6

arguments:
  - -t
  - $(runtime.cores)

stdout: "$(inputs.target.nameroot)_$((Array.isArray(inputs.query) ? inputs.query[0] : inputs.query).nameroot).paf"

outputs:
  alignments: stdout

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-https.rdf
  - https://edamontology.org/EDAM_1.20.owl
