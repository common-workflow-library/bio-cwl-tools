#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/nanoplot:1.29.0--py_0
  SoftwareRequirement:
    packages:
      sra-tools:
        specs: [ "https://github.com/wdecoster/NanoPlot/releases" ]
        version: [ "1.29.0" ]

baseCommand: NanoPlot

inputs:
  max_length:
    type: int?
    inputBinding:
      prefix: '--maxlength'
  min_length:
    type: int?
    inputBinding:
      prefix: '--minlength'
  min_quality:
    type: int?
    inputBinding:
      prefix: '--minqual'
  drop_outliers:
    type: boolean?
    inputBinding:
      prefix: '--drop_outliers'
  log_length:
    type: boolean?
    inputBinding:
      prefix: '--loglength'
  percent_quality:
    type: boolean?
    inputBinding:
      prefix: '--percentqual'
  aligned_length:
    type: boolean?
    inputBinding:
      prefix: '--alength'
  barcoded:
    type: boolean?
    inputBinding:
      prefix: '--barcoded'
  downsample:
    type: int?
    inputBinding:
      prefix: '--downsample'
  run_until:
    type: int?
    inputBinding:
      prefix: '--runtime_until'
  read_type:
    type:
      - "null"
      - type: enum
        symbols: [1D, 2D, 1D2]
    inputBinding:
      prefix: '--readtype'

  color:
    type: string?
    inputBinding:
      prefix: '--color'
  colormap:
    type: string?
    inputBinding:
      prefix: '--colormap'
  format:
    type:
      - type: enum
        symbols: [eps,jpeg,jpg,pdf,pgf,png,ps,raw,rgba,svg,svgz,tif,tiff]
      - "null"
    inputBinding:
      prefix: '--format'
  plots:
    type:
      - type: array
        items:
          type: enum
          symbols: [kde,hex,dot,pauvre]
      - "null"
    inputBinding:
      prefix: '--plots'
  listcolors:
    type: boolean?
    inputBinding:
      prefix: '--listcolors'
  listcolormaps:
    type: boolean?
    inputBinding:
      prefix: '--listcolormaps'
  hide_n50:
    type: boolean?
    inputBinding:
      prefix: '--no-N50'
  show_n50:
    type: boolean?
    inputBinding:
      prefix: '--N50'
  plot_title:
    type: string?
    inputBinding:
      prefix: '--title'
  font_scale:
    type: float?
    inputBinding:
      prefix: '--font_scale'
  dpi:
    type: int?
    inputBinding:
      prefix: '--dpi'
  hide_stats:
    type: boolean?
    inputBinding:
      prefix: '--hide_stats'

  fastq_files:
    type: File[]?
    inputBinding:
      prefix: '--fastq'
  fasta_files:
    type: File[]?
    format: edam:format_1931  # FASTA
    inputBinding:
      prefix: '--fasta'
  rich_fastq_files:
    type: File[]?
    format: edam:format_1930  # FASTQ
    inputBinding:
      prefix: '--fastq_rich'
  minimal_fastq_files:
    type: File[]?
    format: edam:format_1930  # FASTQ
    inputBinding:
      prefix: '--fastq_minimal'
  summary_files:
    type: File[]?
    inputBinding:
      prefix: '--summary'
  bam_files:
    type: File[]?
    format: edam:format_2572  # BAM
    inputBinding:
      prefix: '--bam'
  ubam_files:
    type: File[]?
    inputBinding:
      prefix: '--ubam'
  cram_files:
    type: File[]?
    format: edam:format_3462  # CRAM
    inputBinding:
      prefix: '--cram'
  use_pickle_file:
    type: boolean?
    inputBinding:
      prefix: '--pickle'

outputs:
  dynamic_histogram_read_length:
    type: File
    outputBinding:
      glob: Dynamic_Histogram_Read_length.html
  histogram_read_length:
    type: File
    outputBinding:
      glob: HistogramReadlength.*
  length_v_qual_scatter_plot_dot:
    type: File
    outputBinding:
      glob: LengthvsQualityScatterPlot_dot.*
  length_v_qual_scatter_plot_kde:
    type: File
    outputBinding:
      glob: LengthvsQualityScatterPlot_kde.*
  log_transformed_histogram_read_length:
    type: File
    outputBinding:
      glob: LogTransformed_HistogramReadlength.*
  report:
    type: File
    outputBinding:
      glob: NanoPlot-report.html
  logfile:
    type: File
    outputBinding:
      glob: NanoPlot_*.log
  nanostats:
    type: File
    outputBinding:
      glob: NanoStats.txt
  weighted_histogram_read_length:
    type: File
    outputBinding:
      glob: Weighted_HistogramReadlength.*
  weighted_log_transform_histogram_read_length:
    type: File
    outputBinding:
      glob: Weighted_LogTransformed_HistogramReadlength.*
  yield_by_length_img:
    type: File
    outputBinding:
      glob: Yield_By_Length.*

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
