#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/hopach:v0.0.6
- class:  SoftwareRequirement
  packages:
    hopach:
      specs: [ "http://identifiers.org/biotools/hopach" ]
      version: [ "0.0.6" ]

inputs:

  expression_files:
    type: File[]
    inputBinding:
      prefix: "--input"
    doc: "Input CSV/TSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns"

  expression_aliases:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--name"
    doc: "Input aliases, the order corresponds to --input order. Default: basename of --input files"

  target_column:
    type: string?
    inputBinding:
      prefix: "--target"
    doc: "Target column to be used by hopach clustering. Default: Rpkm"

  combine:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--combine"
    doc: "Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand"

  cluster_method:
    type:
      - "null"
      - type: enum
        symbols: ["row", "column", "both"]
    inputBinding:
      prefix: "--method"
    doc: "Cluster method. Default: both"

  row_dist_metric:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor"]
    inputBinding:
      prefix: "--rowdist"
    doc: "Distance metric for row clustering. Default: cosangle"

  col_dist_metric:
    type:
      - "null"
      - type: enum
        symbols: ["cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor"]
    inputBinding:
      prefix: "--coldist"
    doc: "Distance metric for column clustering. Default: euclid"

  row_logtransform:
    type: boolean?
    inputBinding:
      prefix: "--rowlogtransform"
    doc: "Log2 transform input data prior to running row clustering. Default: false"

  col_logtransform:
    type: boolean?
    inputBinding:
      prefix: "--collogtransform"
    doc: "Log2 transform input data prior to running column clustering. Default: false"

  row_center:
    type:
      - "null"
      - type: enum
        symbols: ["mean", "median"]
    inputBinding:
      prefix: "--rowcenter"
    doc: "Center rows prior to running row clustering. Default: not centered"

  col_center:
    type:
      - "null"
      - type: enum
        symbols: ["mean", "median"]
    inputBinding:
      prefix: "--colcenter"
    doc: "Center columns prior to running column clustering. Default: not centered"

  row_normalize:
    type: boolean?
    inputBinding:
      prefix: "--rownorm"
    doc: "Normalize rows prior to running row clustering. Default: not normalized"

  col_normalize:
    type: boolean?
    inputBinding:
      prefix: "--colnorm"
    doc: "Normalize columns prior to running column clustering. Default: not normalized"

  row_min:
    type: float?
    inputBinding:
      prefix: "--rowmin"
    doc: "Exclude rows from clustering by the min value of a target column. Default: 0"

  keep_discarded:
    type: boolean?
    inputBinding:
      prefix: "--rowkeep"
    doc: "Append excluded rows to the output table after clustering is finished. Default: false"    

  palette:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--palette"
    doc: "Palette color names. Default: red, black, green"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: hopach"

outputs:

  clustering_results:
    type: File
    outputBinding:
      glob: "*_clustering.tsv"
    doc: "Hopach clustering results"

  heatmap_png:
    type: File
    outputBinding:
      glob: "*_heatmap.png"
    doc: "Heatmap ordered by hopach clustering results"

  column_clustering_labels:
    type: File?
    outputBinding:
      glob: "*_column_clustering_labels.tsv"
    doc: "Hopach column clustering labels"

  row_distance_matrix_png:
    type: File?
    outputBinding:
      glob: "*_row_dist_matrix.png"
    doc: "Row distance matrix"

  col_distance_matrix_png:
    type: File?
    outputBinding:
      glob: "*_column_dist_matrix.png"
    doc: "Column distance matrix"

  stderr_log:
    type: File
    outputBinding:
      glob: "hopach_stderr.log"
    doc: "Hopach stderr log"

  stdout_log:
    type: File
    outputBinding:
      glob: "hopach_stdout.log"
    doc: "Hopach stdout log"

baseCommand: ["hopach_order.R"]
stderr: hopach_stderr.log
stdout: hopach_stdout.log

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "hopach"
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Runs hopach clustering algorithm with the combined by specific columns input files.
  Works with minimum two genelist files. The HOPACH clustering algorithm builds a
  hierarchical tree of clusters by recursively partitioning a data set,while ordering
  and possibly collapsing clusters at each level.