#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa:0.7.17--ha92aebf_3
  SoftwareRequirement:
    packages:
      bwa:
        version: [ "0.7.17" ]
        specs: [ https://identifiers.org/biotools/bwa ]

inputs:
  index:
    type: File
    secondaryFiles:
      - ^.amb
      - ^.ann
      - ^.pac
      - ^.sa
    inputBinding:
      position: 200
      valueFrom: $(self.dirname)/$(self.nameroot)
    
  input_files:
    type: File[]
    format:
      - edam:format_1930 # FASTQ (no quality score encoding specified)
      - edam:format_1931 # FASTQ-Illumina
      - edam:format_1932 # FASTQ-Sanger
      - edam:format_1933 # FASTQ-Solexa
    inputBinding:
      position: 201
    
#Optional arguments

  threads:
    type: int?
    inputBinding:
      prefix: "-t"

  min_seed_len:
    type: int?
    inputBinding:
      prefix: "-k"
  
  band_width:
    type: int?
    inputBinding:
      prefix: "-w"

  z_dropoff:
    type: int?
    inputBinding:
      prefix: "-d"

  seed_split_ratio:
    type: float?
    inputBinding:
      prefix: "-r"
    
  max_occ:
    type: int?
    inputBinding:
      prefix: "-c"

  match_score:
    type: int?
    inputBinding:
      prefix: "-A"

  mm_penalty:
    type: int?
    inputBinding:
      prefix: "-B"

  gap_open_penalty:
    type: int?
    inputBinding:
      prefix: "-O"

  gap_ext_penalty:
    type: int?
    inputBinding:
      prefix: "-E"

  clip_penalty:
    type: int?
    inputBinding:
      prefix: "-L"

  unpair_pen:
    type: int?
    inputBinding:
      prefix: "-U"

  rg_line:
    type: string?
    inputBinding:
      prefix: "-R"

  verbose_level:
    type: int?
    inputBinding:
      prefix: "-v"

  is_out_sec_align:
    type: boolean?
    inputBinding:
      prefix: "-a"

  is_mark_short_split:
    type: boolean?
    inputBinding:
      prefix: "-M"

  is_use_hard_clip:
    type: boolean?
    inputBinding:
      prefix: "-H"

  is_multiplexed_pair:
    type: boolean?
    inputBinding:
      prefix: "-p"
      

baseCommand: [bwa, mem]
    
stdout: unsorted_reads.sam

outputs:
  reads_stdout:
    type: stdout
    format: edam:format_2573 # SAM format
    
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
