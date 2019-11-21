#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Aligns reads from ATAC-seq or ChIP-seq to an indexed reference genome

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 30000
  DockerRequirement:
    dockerPull: kerstenbreuer/bowtie2:2.2.6-2
  SoftwareRequirement:
    packages:
      bowtie2:
        specs: [ "http://identifiers.org/biotools/bowtie2" ]
        version: [ "2.3.0" ]

baseCommand: ["bowtie2"]
arguments:
  - valueFrom: --very-sensitive
    position: 1
  - valueFrom: $(runtime.cores) # set the number of threads
    prefix: "-p"
    position: 1
  - position: 10 # prefix for fastq1, differs for paired/single end
    valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "-1";
        }
        else {
          return "-U";
        }
      }
  - valueFrom: $(inputs.fastq1.nameroot).sam
    prefix: "-S"
    position: 6
stderr: $(inputs.fastq1.nameroot).bowtie2_stderr # log file
  

inputs:
  reference_index:
    doc: path to the FM-index files for the chosen reference genome
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
    inputBinding:
      position: 2
      prefix: "-x"
      valueFrom: $(self.path.replace(/\.fa/i,""))
  fastq1:
    type: File
    inputBinding:
      position: 11
  is_paired_end:
    type: boolean
    default: True
  fastq2:
    type: File?
    inputBinding:
      valueFrom: |
        ${
            if ( inputs.is_paired_end ){
                return self;
            }
            else {
              return null;
            }
        }  
      position: 12
      prefix: "-2"
  max_mapping_insert_length:
    doc: usefull for very long fragments, as expected for ATAC
    type: long?
    default: 2000
    inputBinding:
      prefix: --maxins
      position: 1

      
outputs:
  sam:
    type: File
    outputBinding:
      glob: "*.sam"
  bowtie2_log:
    type: stderr
    
