#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: {$include: docker_container.txt}
  SoftwareRequirement:
    packages:
      {$import: software_requirement.yml}

doc: trim primers from mapped reads

baseCommand: ["ivar", "trim"]

inputs:
  bam:
    doc: aligned and sorted reads to be trimmed in BAM format
    type: File
    format: edam:format_2572
    secondaryFiles: .bai
    inputBinding:
      prefix: -i
      position: 10
  bed:
    doc: BED file with primer sequences and positions
    type: File
    format: edam:format_3003
    inputBinding:
      prefix: -b
      position: 20
  output_prefix:
    doc: Prefix for the output BAM file
    type: string
    inputBinding:
      prefix: -p
      position: 30
  min_length:
    doc: Minimum length of read to retain after trimming
    type: int?
    default: 30
    inputBinding:
      position: 1
      prefix: -m
  min_quality:
    doc: Minimum quality threshold for sliding window to pass
    type: int?
    default: 20
    inputBinding:
      position: 1
      prefix: -q
  window_size:
    doc: Width of sliding window
    type: int?
    default: 4
    inputBinding:
      position: 1
      prefix: -s
  primerless_reads:
    doc: Include reads with no primers. By default, reads with no primers are excluded
    type: boolean
    default: false
    inputBinding:
      prefix: -e

outputs:
  bam_filtered:
    type: File
    format: edam:format_2572
    outputBinding:
      glob: $(inputs.output_prefix).bam
  
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
