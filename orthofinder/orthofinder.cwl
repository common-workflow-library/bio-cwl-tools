#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: "OrthoFinder"
doc: |
  
  OrthoFinder: phylogenetic orthology inference for comparative genomics

  SIMPLE USAGE:
  Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

  Source: https://github.com/davidemms/OrthoFinder

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: 
      - $(inputs.sequences)

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/orthofinder:2.5.5--hdfd78af_2
  SoftwareRequirement:
    packages:
      - package: OrthoFinder
        version: [ "2.5.5" ]
        specs: [ https://identifiers.org/rrid/RRID:SCR_017118 ]    
  ResourceRequirement:
    ramMin: $(16 * 1024)
    coresMin: 4

inputs:
  sequences:
    type: File[]
    format: edam:format_1929
    label: "FASTA sequence files"
    doc: |
      "OrthoFinder requires a writeable directory of sequence files as input (`-f <dir>`). 
      Hence the files are staged via `InitialWorkDirRequirement` and passed as $(runtime.outdir) to argument -f."

outputs:
  OrthoFinderOutDir:
    type: Directory
    doc: "OrthoFinder outputs a directory OrthoFinder/Results_<date in MMMDD format>"
    outputBinding:
      glob: OrthoFinder/Results*

baseCommand: orthofinder

arguments:
  - prefix: -f
    valueFrom: $(runtime.outdir)

$namespaces:
  edam: https://edamontology.org/
  s: https://schema.org/
$schemas:
  - https://edamontology.org/EDAM_1.25.owl
  - https://schema.org/version/latest/schemaorg-current-https.rdf
