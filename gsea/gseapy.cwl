#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/gseapy:0.13.0--py310hbee2dd9_0
  SoftwareRequirement:
    packages:
      gseapy:
        specs:
          - https://anaconda.org/bioconda/gseapy
          - https://github.com/zqfang/GSEApy

inputs:
  read_counts_file:
    type: File
    format:
      - iana:text/plain
      - edam:format_3709  # GCT
    inputBinding:
      prefix: "-d"
    doc: "Input gene expression dataset file in txt or gct format. Same with GSEA"

  phenotypes_file:
    type: File
    inputBinding:
      prefix: "-c"
    doc: "Input class vector (phenotype) file in CLS format. Same with GSEA"

  gene_set_database:
    type:
    - File
    - type: enum
      name: "genesetdatabase"
      symbols:
      - H_hallmark_gene_sets
      - C1_positional_gene_sets
      - C2_curated_gene_sets
      - C3_regulatory_target_gene_sets
      - C4_computational_gene_sets
      - C5_GO_gene_sets
      - C6_oncogenic_signatures
      - C7_immunologic_signatures
    inputBinding:
      prefix: "-g"
    doc: "Gene set database"

  permutation_type:
    type:
    - "null"
    - type: enum
      name: "permutationtype"
      symbols:
      - gene_set
      - phenotype
    inputBinding:
      prefix: "-t"
    doc: "Permutation type. Default: gene_set"

  permutation_count:
    type: int?
    inputBinding:
      prefix: "-n"
    doc: "Number of random permutations. For calculating esnulls. Default: 1000"

  min_gene_set_size:
    type: int?
    inputBinding:
      prefix: "--min-size"
    doc: "Min size of input genes presented in Gene Sets. Default: 15"

  max_gene_set_size:
    type: int?
    inputBinding:
      prefix: "--max-size"
    doc: "Max size of input genes presented in Gene Sets. Default: 500"

  ranking_metrics:
    type:
    - "null"
    - type: enum
      name: "rankingmetrics"
      symbols:
      - signal_to_noise
      - t_test
      - ratio_of_classes
      - diff_of_classes
      - log2_ratio_of_classes
    inputBinding:
      prefix: "-m"
    doc: "Methods to calculate correlations of ranking metrics. Default: log2_ratio_of_classes"

  ascending_rank_sorting:
    type: boolean?
    inputBinding:
      prefix: "-a"
    doc: "Ascending rank metric sorting order. Default: False"

  graphs_count:
    type: int?
    inputBinding:
      prefix: "--graph"
    doc: "Numbers of top graphs produced. Default: 20"

  seed:
    type: int?
    inputBinding:
      prefix: "-s"
    doc: "Number of random seed. Default: None"

  threads:
    type: int?
    inputBinding:
      prefix: "-p"
    doc: "Threads number"

outputs:
  enrichment_report:
    type: File?
    format: iana:text/csv
    outputBinding:
      glob: "GSEApy_reports/*.report.csv"

  enrichment_plots:
    type: File[]
    format: iana:application/pdf
    outputBinding:
      glob: "GSEApy_reports/*.gsea.pdf"

  enrichment_heatmaps:
    type: File[]
    format: iana:application/pdf
    outputBinding:
      glob: "GSEApy_reports/*.heatmap.pdf"

baseCommand: [gseapy, gsea]

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

label: GSEApy
s:license: http://www.apache.org/licenses/LICENSE-2.0

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
  GSEAPY: Gene Set Enrichment Analysis in Python
  ==============================================

  Gene Set Enrichment Analysis is a computational method that determines whether an a priori
  defined set of genes shows statistically significant, concordant differences between two
  biological states (e.g. phenotypes).

  GSEA requires as input an expression dataset, which contains expression profiles for multiple samples.
  While the software supports multiple input file formats for these datasets, the tab-delimited GCT format
  is the most common. The first column of the GCT file contains feature identifiers (gene ids or symbols in
  the case of data derived from RNA-Seq experiments). The second column contains a description of the feature;
  this column is ignored by GSEA and may be filled with “NA”s. Subsequent columns contain the expression
  values for each feature, with one sample's expression value per column. It is important to note that there
  are no hard and fast rules regarding how a GCT file's expression values are derived. The important point is
  that they are comparable to one another across features within a sample and comparable to one another
  across samples. Tools such as DESeq2 can be made to produce properly normalized data (normalized counts)
  which are compatible with GSEA.
