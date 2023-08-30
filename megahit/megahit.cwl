#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  forward_reads:
    type: File
    format: edam:format_1930  # FASTQ
  reverse_reads:
    type: File
    format: edam:format_1930  # FASTQ

baseCommand: megahit

arguments:
 - -t
 - $(runtime.cores)
 - "-1"
 - $(inputs.forward_reads.path)
 - "-2"
 - $(inputs.reverse_reads.path)
 - -o
 - $(runtime.outdir)/results

outputs:
  metagenome_assembly:
    type: File
    format: edam:format_1929  # FASTA
    outputBinding:
     glob: $(runtime.outdir)/results/final.contigs.fa

hints:
  SoftwareRequirement:
    packages:
      megahit:
        specs:
          - https://identifiers.org/biotools/megahit
          - https://identifiers.org/rrid/RRID:SCR_018551
  DockerRequirement:
    dockerPull: quay.io/biocontainers/megahit:1.2.9--h43eeafb_4
  ResourceRequirement:
    coresMin: 2

$namespaces:
  edam: https://edamontology.org
