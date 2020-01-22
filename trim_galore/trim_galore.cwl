#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Adaptor trimming of reads (single or paired end) in fastq format.

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 7000
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7

baseCommand: trim_galore

inputs:
  # main input
  fastq1:
    doc: |
      raw reads in fastq format; can be gzipped;
      if paired end, the file contains the first reads;
      if single end, the file contains all reads
    type: File
    inputBinding:
      position: 10
  fastq2:
    doc: |
      (optional) raw reads in fastq format; can be gzipped;
      if paired end, the file contains the second reads;
      if single end, the file does not exist
    type: File?
    inputBinding:
      position: 11
  adapter1:
    doc: |
      Adapter sequence for first reads.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the string to "illumina", "nextera", or "small_rna".
    type: string?
  adapter2:
    doc: |
      Adapter sequence for second reads - only for paired end data.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the adapter1 string to "illumina", "nextera", or "small_rna".
    type: string?
    
  # additional optional input:
  qual_trim_cutoff:
    doc: trim all base with a phred score lower than this valueFrom
    type: int
    default: 20
    inputBinding:
      prefix: --quality
      position: 1
  min_read_length:
    doc: discard reads that get shorter than this value
    type: int
    default: 20
    inputBinding:
      prefix: --length
      position: 1
  min_unpaired_read_rescue_length:
    doc: |
      if only one read of a pair passes the qc and adapter trimming,
      it needs at least this length to be rescued
    type: int
    default: 35
  min_adapter_overlap:
    doc: minimum overlap with adapter seq in bp needed to trim
    type: int
    default: 1
    inputBinding:
      prefix: --stringency 
      position: 1

arguments:

  ## hard-coded parameters:
  - prefix: --fastqc_args
    valueFrom: "\"--noextract\""
    position: 1
    # fastqc data remains zipped
  - prefix: --gzip
    position: 1
    # gzip output fastq
  

  ## variable arguments:
  - valueFrom: |
      ${
        if ( inputs.adapter1 == "illumina" ){ return "--illumina" }
        else if ( inputs.adapter1 == "nextera" ){ return "--nextera" }
        else if ( inputs.adapter1 == "small_rna" ){ return "--small_rna" }
        else { return null }
      }
    position: 1
  - prefix: --adapter
    valueFrom: |
      ${
        if ( inputs.apdater1 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter1
        } else {
          return null
        }
      }
    position: 1
  - prefix: --adapter2
    valueFrom: |
      ${
        if ( inputs.fastq2 != null && inputs.apdater2 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter2
        } else {
          return null
        }
      }
    position: 1
  - valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return "--paired" }
      }
    position: 1
    # turn on paired mode
  - valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return "--retain_unpaired" }
      }
    position: 1
    # if only one read of a read pair becomes to short,
    # it will be written to .unpaired_1.fq or .unpaired_2.fq
  - prefix: --length_1
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return inputs.min_unpaired_read_rescue_length }
      }
    position: 1
    # min bp length for a read to be written to .unpaired_1.fq
  - prefix: --length_2
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return inputs.min_unpaired_read_rescue_length }
      }
    position: 1
    # min bp length for a read to be written to .unpaired_2.fq

outputs:
  fastq1_trimmed:
    type: File
    outputBinding:
      glob: |
        ${
            if ( inputs.fastq2 == null  ){ return "*trimmed.fq*" }
            else { return "*val_1.fq*" }
        }
  fastq2_trimmed:    
    type: File?
    outputBinding:
      glob: "*val_2.fq*"
  fastq1_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: "*unpaired_1.fq*"
  fastq2_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: "*unpaired_2.fq*"
  trim_galore_log: # can be used by multiqc
    type:
      type: array # since one or two matches (single/paired end)
      items: File
    outputBinding:
      glob:  "*trimming_report.txt"
  trimmed_fastqc_html:
    doc: html report of post-trimming fastqc
    type:
      type: array # since one or two matches (single/paired end)
      items: File
    outputBinding:
      glob: "*fastqc.html"
  trimmed_fastqc_zip:
    doc: all data of post-trimming fastqc e.g. figures
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputBinding:
      glob: "*fastqc.zip"
