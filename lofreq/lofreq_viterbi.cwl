#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: "viterbi: Viterbi realignment"
doc: |
  Probabilistic realignment of your already mapped reads, which corrects
  mapping errors (run right after mapping). Not recommended for non-Illumina
  data.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/lofreq:2.1.4--py27hc3dfafe_1

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference)
    
baseCommand: [lofreq, viterbi]

arguments:
  - prefix: --out
    valueFrom: $(inputs.reads.nameroot)_realigned.bam

inputs:
  reference:
    type: File
    format: edam:format_1930  # FASTA
    inputBinding:
      prefix: --ref
      
  reads:  
    type: File
    format: edam:format_2572  # BAM
    inputBinding: {}

  keepflags:
    type: boolean?
    label: Don't delete flags MC, MD, NM, and A?
    doc: |
      These flags are all prone to getting invalidated during realignment.
      Keep them only if you know what you are doing.
    inputBinding:
      prefix: --keepflags
    default: false
      
  defqual:
    type: int?
    inputBinding:
      prefix: --defqual

outputs:
  realigned:
    type: File
    format: edam:format_2572  # BAM
    outputBinding:
      glob: $(inputs.reads.nameroot)_realigned.bam

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
