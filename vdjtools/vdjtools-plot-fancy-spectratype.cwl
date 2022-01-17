class: CommandLineTool
cwlVersion: v1.1


requirements:
- class: InlineJavascriptRequirement
- class: ResourceRequirement
  ramMin: 3814
  coresMin: 2
- class: DockerRequirement
  dockerPull: yyasumizu/vdjtools
  

inputs:

  top:
    type: int?
    inputBinding:
      prefix: "--top"
      position: 5
    doc: |
      Number of top clonotypes to visualize. Should not exceed 20.
      Default: 10

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
  
  fancy_spectratype_file:
    type: File
    outputBinding:
      glob: "*.txt"

  fancy_spectratype_plot:
    type: File
    outputBinding:
      glob: "*.pdf"


baseCommand: ["vdjtools", "PlotFancySpectratype"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "VDJtools Plot Spectratype"
s:name: "VDJtools Plot Spectratype"
s:alternateName: "Plots a spectratype that displays CDR3 lengths for top N clonotypes in a given sample"

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
  VDJtools is an open-source Java/Groovy-based framework designed to facilitate analysis of
  immune repertoire sequencing (RepSeq) data. VDJtools computes a wide set of statistics and
  is able to perform various forms of cross-sample analysis. Both comprehensive tabular output
  and publication-ready plots are provided.
  
  Plots a spectratype that also displays CDR3 lengths for top N clonotypes in a given sample.
  This plot allows to detect the highly-expanded clonotypes.