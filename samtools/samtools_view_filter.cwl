#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

doc: for single end data

baseCommand: ["samtools", "view"]
arguments:
  - valueFrom: -c
    position: 1
    # count reads
  - valueFrom: "4"
    prefix: -F
    position: 1
    # exclude unmapped
  # when paired end, "-f 3" is added to the command line
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "-f";
        }
        else {
          return null;
        }
      }
    position: 2
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "3";
        }
        else {
          return null;
        }
      }
    position: 3
stdout: $(inputs.bam.nameroot)_filt.bam

inputs:
  bam:
    doc: aligned reads to be checked in bam format
    type: File
    inputBinding:
      position: 10
  is_paired_end:
    doc: if paired end, only properly paired reads pass
    type: boolean
    default: true
  min_mapping_quality:
    doc: Reads with a mapping quality below this will be excluded
    type: int
    default: 20
    inputBinding:
      position: 1
      prefix: -q

outputs:
  bam_filtered:
    type: stdout
  
  
