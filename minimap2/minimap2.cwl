#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: minimap2

doc: This CWL file defines running minimap2 to align some sequences to a database.
  We assume the database has been indexed. This is not necessary but we will do it
  in our use case

hints:
  DockerRequirement:
    dockerPull: ttubb/minimap2:release-0.2.0
  SoftwareRequirement:
    packages:
      minimap2:
        specs: [ https://github.com/lh3/minimap2 ]
        version: [ "2.16" ]

inputs:
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: "-t"
  preset:
    type: string
    default: sr
    inputBinding:
      position: 2
      prefix: "-x"
  samOutput:
    type: boolean
    default: true
    inputBinding:
      position: 3
      prefix: "-a"
  indexFile:
    type: File
    inputBinding:
      position: 5
  fastqFiles:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      position: 6

arguments:
  - -o
  - output.sam

outputs:
  samfile:
    type: File
    outputBinding:
      glob: "output.sam"
    streamable: true

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.20.owl
