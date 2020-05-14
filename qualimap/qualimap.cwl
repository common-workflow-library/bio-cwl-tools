#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: qualimap-qc
doc: |-
  This is qualimap CWL tool definition http://qualimap.bioinfo.cipf.es/.
  It perform RNA-seq QC analysis on paired-end data http://qualimap.bioinfo.cipf.es/doc_html/command_line.html.

hints:
  ResourceRequirement:
    ramMin: 4000
    coresMin: 1
  SoftwareRequirement:
    packages:
      qualimap:
        version: [ "2.2.2d" ]
        specs:
          - https://identifiers.org/rrid/RRID:SCR_001209
          - https://bio.tools/qualimap
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/qualimap:2.2.2d--1'
requirements:
  - class: ShellCommandRequirement

# Base command
baseCommand: [ qualimap, rnaseq ]

arguments:
  - valueFrom: |-
      --paired --java-mem-size="$(inputs.javamem)"
    shellQuote: false
  - prefix: -outdir
    valueFrom: $(inputs.inputBam.nameroot)

inputs:
  javamem:
    type: string
    default: "1200M"
    label: Set desired Java heap memory size

  outdir:
    type: string
    inputBinding:
      prefix: "-outdir"
    doc: |
      Output folder for HTML report and raw data.

  algo:
    type:
      - "null"
      - type: enum
        symbols:
          - uniquely-mapped-reads
          - proportional
    inputBinding:
      prefix: "--algorithm"
    label: Counting algorithm

  inputBam:
    type: File
    secondaryFiles: .bai
    inputBinding:
      prefix: "-bam"
    doc: |
      Input mapping file in BAM format.

  seqProtocol:
    type:
      - "null"
      - type: enum
        symbols:
           - strand-specific-forward
           - strand-specific-reverse
           - non-strand-specific
    inputBinding:
      prefix: "--sequencing-protocol"
    label: Sequencing library protocol
  
  gtf:
    type: File
    inputBinding:
      prefix: "-gtf"
    doc: |
      Region file in GTF, GFF or BED format. 
      If GTF format is provided, counting is based on
      attributes, otherwise based on feature name.

outputs:
  qualimapQC:
    type: Directory
    outputBinding:
      glob: "$(inputs.outdir)" 

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:codeRepository: https://github.com/common-workflow-library/bio-cwl-tools
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:creator:
- class: s:Organization
  s:legalName: "University of Melbourne Centre for Cancer Research"
  s:member:
  - class: s:Person
    s:name: Dr. Sehrish Kanwal
    s:email: mailto:kanwals@unimelb.edu.au
  
