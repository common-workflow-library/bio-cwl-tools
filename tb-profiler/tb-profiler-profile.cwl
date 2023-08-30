#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Profile M. tuberculosis whole genome sequencing samples for drug resistance and lineage
label: TB-Profiler profile

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-6553-5274
    s:email: mailto:pvh@sanbi.ac.za
    s:name: Peter van Heusden

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/tb-profiler:3.0.0--pypyh3252c3a_0
  SoftwareRequirement:
    packages:
      tb-profiler:
        specs: [ https://identifiers.org/biotools/tb-profiler ]
        version: [ "3.0.0" ]

hints:
  ResourceRequirement:
    coresMin: 1
    # TODO: estimate ramMin
  SoftwareRequirement:
    packages:
      tb-profiler:
        version: [ "3.0.0" ]
        specs:
          - https://anaconda.org/bioconda/tb-profiler
          - https://github.com/jodyphelan/TBProfiler/

inputs:
  sequences:
    label: "Input sequences"
    type:
      - type: record
        fields:
          read1:
            label: "First read file"
            type: File
            inputBinding:
              prefix: --read1
              position: 10
          read2:
            label: "Second read file"
            type: File?
            inputBinding:
              prefix: --read2
              position: 20
      - type: record
        fields:
          bam:
            label: "Aligned reads in BAM format"
            doc: "Make sure it has been generated using the H37RV reference genome (GCA_000195955.2)"
            type: File
            inputBinding:
              prefix: --bam
              position: 10
  database:
    label: "Mutation database to use for drug resistance prediction"
    type:
      - type: record
        fields:
          panel_name:
            label: "Mutation panel name"
            type: string
            inputBinding:
              prefix: --db
              position: 30
      - type: record
        fields:
          external_db:
            label: "Path prefix to db files"
            type: string
            inputBinding:
              prefix: --external_db
              position: 30
    default: 
      panel_name: "tbdb"
  platform:
    label: "NGS platform used to generate data"
    type:
      type: enum
      symbols:
        - illumina
        - nanopore
      inputBinding:
        prefix: --platform
        position: 40      
    default: "illumina"
  output_prefix:
    label: "Prefix to use for all results generated"
    type: string?
    default: "tbprofiler"
    inputBinding:
      prefix: --prefix
      position: 50
  csv_output:
    label: "Produce CSV output"
    type: boolean?
    inputBinding:
      prefix: --csv
      position: 50
  text_output:
    label: "Produce text output"
    type: boolean?
    inputBinding:
      prefix: --txt
      position: 50
  pdf_output:
    label: "Produce PDF output (using pdflatex)"
    type: boolean?
    inputBinding:
      prefix: --pdf
      position: 50
  mapper:
    label: "Mapping tool to use."
    doc: "Defaults to bwa if the illumina sequencing platform is selected, and defaults to minimap2 for nanopore"
    type:
      - type: enum
        symbols:
          - "bwa"
          - "minimap2"
          - "bowtie2"
          - "bwa-mem2"
        inputBinding:
          prefix: --mapper
          position: 60
      - "null"      
  caller:
    label: "Variant calling tool to use"
    doc: "Default is bcftools"
    type:
      - type: enum
        symbols:
          - "bcftools"
          - "gatk"
          - "freebayes"
        inputBinding:
          prefix: --caller
          position: 60
      - "null"
  calling_parameters:
    label: "Override default parameters for variant calling"
    type: string?
    inputBinding:
      prefix: --calling_params
      position: 60
  min_depth:
    label: "Minimum depth required to call variant"
    doc: "Bases with depth below this cutoff will be marked as missing"
    type: int?
    inputBinding:
      prefix: --min_depth
      position: 60
  allele_frequency:
    label: "Minimum allele frequency to call variants"
    type: float?
    inputBinding:
      prefix: --af
      position: 60
  dr_allele_frequency:
    label: "Minimum allele frequency to use variants for prediction"
    type: float?
    inputBinding:
      prefix: --reporting_af
      position: 60
  threads:
    label: "Threads to use"
    type: int?
    inputBinding:
      prefix: --threads
      position: 60
  no_trim:
    label: "Don't trim reads using trimmomatic"
    type: boolean?
    inputBinding:
      prefix: --no_trim
      position: 60
  no_flagstat:
    label: "Don't collect flagstats"
    type: boolean?
    inputBinding:
      prefix: --no_flagstat
      position: 60
  no_delly:
    label: "Don't run delly"
    type: boolean?
    inputBinding:
      prefix: --no_delly
      position: 60
  
  coverage_fraction_threshold:
    label: "Cutoff used to calculate fraction of region covered by <= this value"
    type: float?
    inputBinding:
      prefix: --coverage_fraction_threshold
      position: 60

outputs:
  json_output_file:
    type: File
    outputBinding:
      glob: results/$(inputs.output_prefix)*.json
  csv_output_file:
    type: File?
    outputBinding:
      glob: results/$(inputs.output_prefix)*.csv
  text_output_file:
    type: File?
    outputBinding:
      glob: results/$(inputs.output_prefix)*.txt
  pdf_output_file:
    type: File?
    outputBinding:
      glob: results/$(inputs.output_prefix)*.pdf


baseCommand: [ "tb-profiler", "profile" ]

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-https.rdf
