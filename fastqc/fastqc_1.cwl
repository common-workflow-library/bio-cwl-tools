#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Run fastqc on raw reads in FASTQ format (single or paired end) or aligned reads in BAM.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 5000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1
  SoftwareRequirement:
    packages:
      fastqc:
        specs: [ https://identifiers.org/biotools/fastqc ]
        version: [ "0.11.9" ]


baseCommand: "fastqc"
arguments: 
  - valueFrom: $(runtime.outdir)
    prefix: "-o"
  - valueFrom: "--noextract"

inputs:
  fastq1:
    type: File?
    inputBinding:
      position: 1
  fastq2:
    type: File?
    inputBinding:
      position: 2
  bam:
    type: File?
    inputBinding:
      position: 1
 
outputs:
  fastqc_zip:
    doc: all data e.g. figures
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.zip"
  fastqc_html:
    doc: html report showing results from zip
    type:
      type: array
      items: File
    outputBinding:
      glob: "*_fastqc.html"
    
