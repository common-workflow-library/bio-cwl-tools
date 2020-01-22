#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.8
- class: SoftwareRequirement
  packages:
    deseq:
      specs: [ "http://identifiers.org/biotools/deseq" ]
      version: [ "0.0.8" ] ## TODO: Update!

inputs:

  untreated_files:
    type:
      - File
      - File[]
    inputBinding:
      position: 5
      prefix: "-u"
    doc: |
      Untreated input CSV/TSV files

  treated_files:
    type:
      - File
      - File[]
    inputBinding:
      position: 6
      prefix: "-t"
    doc: |
      Treated input CSV/TSV files

  untreated_col_suffix:
    type: string?
    inputBinding:
      position: 7
      prefix: "-un"
    doc: |
      Suffix for untreated RPKM column name

  treated_col_suffix:
    type: string?
    inputBinding:
      position: 8
      prefix: "-tn"
    doc: |
      Suffix for treated RPKM column name

  output_filename:
    type: string
    inputBinding:
      position: 9
      prefix: "-o"
    doc: |
      Output TSV filename

  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: '-p'
    doc: |
      Run script using multiple threads

outputs:

  diff_expr_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

  plot_lfc_vs_mean:
    type: File
    outputBinding:
      glob: "*001.png"

  gene_expr_heatmap:
    type: File
    outputBinding:
      glob: "*002.png"

baseCommand: [run_deseq.R]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "deseq_advanced"
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
  Tool runs DESeq/DESeq2 script similar to the original one from BioWArdrobe.
  untreated_files and treated_files input files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file in
  the untreated group. Output is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq/DESeq2 always compares untreated_vs_treated groups