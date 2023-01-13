#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

hints:
  ResourceRequirement:
    ramMin: 3814
    coresMin: 2
  DockerRequirement:
    dockerPull: yyasumizu/vdjtools
  
inputs:
  top:
    type: int?
    inputBinding:
      prefix: "--top"
      position: 5

  vdj_file:
    type: File
    inputBinding:
      position: 6

  output_prefix:
    type: string?
    default: "./"
    inputBinding:
      position: 7

outputs:
  quantile_stats_file:
    type: File
    outputBinding:
      glob: "*.txt"

  quantile_stats_plot:
    type: File
    outputBinding:
      glob: "*.pdf"

baseCommand: ["vdjtools", "PlotQuantileStats"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

label: "VDJtools Plot Quantile Stats"
s:alternateName: "Plots a three-layer donut chart to visualize the repertoire clonality"

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
  VDJtools is an open-source Java/Groovy-based framework designed to facilitate analysis of
  immune repertoire sequencing (RepSeq) data. VDJtools computes a wide set of statistics and
  is able to perform various forms of cross-sample analysis. Both comprehensive tabular output
  and publication-ready plots are provided.
  
  Plots a three-layer donut chart to visualize the repertoire clonality.
  - First layer (“set”) includes the frequency of singleton (“1”, met once), doubleton
    (“2”, met twice) and high-order (“3+”, met three or more times) clonotypes. Singleton and
    doubleton frequency is an important factor in estimating the total repertoire diversity,
    e.g. Chao1 diversity estimator (see Colwell et al). We have also recently shown that in
    whole blood samples, singletons have very nice correlation with the number of naive T-cells,
    which are the backbone of immune repertoire diversity.
  - The second layer (“quantile”), displays the abundance of top 20% (“Q1”), next 20% (“Q2”), ...
    (up to “Q5”) clonotypes for clonotypes from “3+” set. In our experience this quantile plot is
    a simple and efficient way to display repertoire clonality.
  - The last layer (“top”) displays the individual abundances of top N clonotypes.
