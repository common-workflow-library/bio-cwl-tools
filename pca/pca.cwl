#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/pca:v0.0.4

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
    doc: "Target column name to be used by PCA. Default: Rpkm"

  combine:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "--combine"
    doc: "Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: pca_"

outputs:

  pca1_vs_pca2_plot:
    type: File
    outputBinding:
      glob: "*001.png"
    doc: "PCA1 vs PCA2 plot"

  pca2_vs_pca3_plot:
    type: File
    outputBinding:
      glob: "*002.png"
    doc: "PCA2 vs PCA3 plot"

  variance_plot:
    type: File
    outputBinding:
      glob: "*003.png"
    doc: "Variance plot"
    
  pca_3d_plot:
    type: File?
    outputBinding:
      glob: "*004.png"
    doc: "First three principal components plot"

  pca_file:
    type: File
    outputBinding:
      glob: "*.tsv"
    doc: "PCA analysis results exported as TSV"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

baseCommand: ["run_pca.R"]
stderr: pca_stderr.log
stdout: pca_stdout.log

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "pca"
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
  Principal Component Analysis
  --------------

  Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert
  a set of observations of possibly correlated variables (entities each of which takes on various numerical values)
  into a set of values of linearly uncorrelated variables called principal components.

  The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix,
  not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy.