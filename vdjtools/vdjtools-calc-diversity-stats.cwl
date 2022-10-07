#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    ramMin: 3814
    coresMin: 2
  DockerRequirement:
    dockerPull: yyasumizu/vdjtools
  
inputs:
  intersect_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "strict"
      - "nt"
      - "ntV"
      - "ntVJ"
      - "aa"
      - "aaV"
      - "aaVJ"
      - "aa!nt"
    inputBinding:
      prefix: "--intersect-type"
      position: 5
    doc: |
      Specifies which clonotype features (CDR3 sequence, V/J segments, hypermutations)
      will be compared when checking if two clonotypes match.
      Default: strict

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
  diversity_stats_file:
    type: File
    outputBinding:
      glob: "*.exact.txt"

baseCommand: ["vdjtools", "CalcDiversityStats"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "VDJtools Calc Diversity Stats"
s:alternateName: "Computes a set of diversity statistics"

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

  Computes a set of diversity statistics, including
  - Observed diversity, the total number of clonotypes in a sample
  - Lower bound total diversity (LBTD) estimates
  - Chao estimate (denoted chao1)
  - Efron-Thisted estimate
  - Diversity indices
    * Shannon-Wiener index. The exponent of clonotype frequency distribution entropy is returned.
    * Normalized Shannon-Wiener index. Normalized (divided by log[number of clonotypes]) entropy
      of clonotype frequency distribution. Note that plain entropy is returned, not its exponent.
    * Inverse Simpson index
  - Extrapolated Chao diversity estimate (denoted chaoE)
  - The d50 index, a recently developed immune diversity estimate

  --downsample-to, --extrapolate-to, and --resample-trials parameters are omitted as the input is
  a single file. Only not resampled output is returned.
