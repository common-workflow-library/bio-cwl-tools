#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Removal of duplicates from aligned reads.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/picard_tools:2.17.4
  
baseCommand: ["java", "-jar"]
arguments:
  - valueFrom: "MarkDuplicates"
    position: 2
  - valueFrom: $(inputs.bam_sorted.nameroot)_duprem.bam
    prefix: "OUTPUT="
    separate: false
    position: 13
  - valueFrom: $(inputs.bam_sorted.nameroot)_duprem.log
    prefix: "METRICS_FILE="
    separate: false
    position: 13
    # log file
  - valueFrom: "REMOVE_DUPLICATES=TRUE"
    position: 14
  - valueFrom: "ASSUME_SORTED=TRUE"
    position: 15
  - valueFrom: "VALIDATION_STRINGENCY=SILENT"
    position: 16
  - valueFrom: "VERBOSITY=INFO"
    position: 17
  - valueFrom: "QUIET=false"
    position: 17
stderr: $(inputs.bam_sorted.nameroot).picard_markdup.log

inputs:
  bam_sorted:
    doc: sorted bam input file
    type: File
    inputBinding:
      prefix: "INPUT="
      separate: false
      position: 11
  path_to_picards:
    type: string
    default: "/bin/picard.jar"
    inputBinding:
      position: 1
      
outputs:
  bam_duprem:
    type: File
    outputBinding:
      glob: $(inputs.bam_sorted.nameroot)_duprem.bam
  picard_markdup_log:
    type: stderr
    
