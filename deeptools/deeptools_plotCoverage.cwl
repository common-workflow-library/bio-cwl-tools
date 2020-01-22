#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Create a plot showing global coverage for reads provided in BAM.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1
  SoftwareRequirement:
    packages:
      deeptools:
        specs: [ "http://identifiers.org/biotools/deeptools" ]
        version: [ "3.1.1" ]

baseCommand: ["plotCoverage"]
arguments:
  - valueFrom: $(inputs.sample_id)
    prefix: --labels
    position: 10
  - valueFrom: $(inputs.sample_id).plot_cov.png
    prefix: --plotFile
    position: 10
  - valueFrom: $(inputs.sample_id).plot_cov.tsv
    prefix: --outRawCounts
    position: 10
  
inputs:
  bam:
    doc: must be indexed
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bamfiles
  extendReads:
    doc: For single end reads, the fragment size should be specified here.
    type: int?
    inputBinding:
      prefix: --extendReads
      position: 10
  sample_id:
    type: string
    
outputs:
  qc_plot_coverage_plot:
    type: File
    outputBinding:
      glob: $(inputs.sample_id).plot_cov.png
  qc_plot_coverage_tsv:
    type: File
    outputBinding:
      glob: $(inputs.sample_id).plot_cov.tsv
    
