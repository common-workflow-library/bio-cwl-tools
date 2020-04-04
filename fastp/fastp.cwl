#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Modified from https://github.com/nigyta/bact_genome/blob/master/cwl/tool/fastp/fastp.cwl
requirements:
    InlineJavascriptRequirement: {}
hints:
    DockerRequirement:
        dockerPull: quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0

baseCommand: [fastp]

arguments:
    - prefix: -o
      valueFrom: $(inputs.fastq1.nameroot).fastp.fastq
    - |
      ${
        if (inputs.fastq2){
          return '-O';
        } else {
          return '';
        }
      }
    - |
      ${
        if (inputs.fastq2){
          return inputs.fastq2.nameroot + ".fastp.fastq";
        } else {
          return '';
        }
      }

inputs:
    fastq1:
      type: File
      format:
        - edam:format_1930 # FASTA
        - edam:format_1931 # FASTQ
      inputBinding:
        prefix: -i
    fastq2:
      format:
        - edam:format_1930 # FASTA
        - edam:format_1931 # FASTQ
      type: File?
      inputBinding:
        prefix: -I
    threads:
      type: int?
      default: 1
      inputBinding:
        prefix: --thread
    qualified_phred_quality:
      type: int?
      default: 20
      inputBinding:
        prefix: --qualified_quality_phred
    unqualified_phred_quality:
      type: int?
      default: 20
      inputBinding:
        prefix: --unqualified_percent_limit
    min_length_required:
      type: int?
      default: 50
      inputBinding:
        prefix: --length_required
    force_polyg_tail_trimming:
      type: boolean?
      inputBinding:
        prefix: --trim_poly_g
    disable_trim_poly_g:
      type: boolean?
      default: true
      inputBinding:
        prefix: --disable_trim_poly_g
    base_correction:
      type: boolean?
      default: true
      inputBinding:
        prefix: --correction
    


outputs:
    out_fastq1:
       type: File
       format: $(inputs.fastq1.format)
       outputBinding:
           glob: $(inputs.fastq1.nameroot).fastp.fastq
    out_fastq2:
       type: File?
       format: $(inputs.fastq2.format)
       outputBinding:
           glob: $(inputs.fastq2.nameroot).fastp.fastq
    html_report:
      type: File
      outputBinding:
        glob: fastp.html
    json_report:
      type: File
      outputBinding:
        glob: fastp.json

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
