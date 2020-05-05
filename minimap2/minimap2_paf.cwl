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
        specs: [ https://github.com/lh3/minimap2 ]
        version: [ "2.17" ]
  ResourceRequirement:
    coresMin: 8
    coresMax: 32
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))

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

stdout: $(inputs.indexFile.nameroot)_$(inputs.fastqFiles[0].nameroot).paf

outputs:
  alignments: stdout

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.20.owl
