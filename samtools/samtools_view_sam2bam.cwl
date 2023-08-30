#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Convert SAM to BAM.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.14--hb421002_0
  SoftwareRequirement:
    packages:
      samtools:
        specs: [ https://identifiers.org/biotools/samtools ]
        version: [ "1.14" ]

baseCommand: ["samtools", "view"]

arguments:
  - valueFrom: -h
    position: 1
    # include the headers
  - valueFrom: -b
    position: 1
    # output in bam format
  - --no-PG  # but don't add our own PG header for more reproduciblity

stdout: $(inputs.sam.nameroot).bam

inputs:
  sam:
    doc: reads to be checked in sam format
    type: File
    format: edam:format_2573  # SAM
    inputBinding:
      position: 2

outputs:
  bam:
    type: File
    outputBinding:
      glob: "*.bam"
    format: edam:format_2572  # BAM

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl 
