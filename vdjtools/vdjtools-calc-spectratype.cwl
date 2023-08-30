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
  unweighted:
    type: boolean?
    inputBinding:
      prefix: "--unweighted"
      position: 5

  amino_acid:
    type: boolean?
    inputBinding:
      prefix: "--amino-acid"
      position: 6

  vdj_file:
    type: File
    inputBinding:
      position: 7

  output_prefix:
    type: string?
    default: "./"
    inputBinding:
      position: 8


outputs:
  
  spectratype_insert_wt_file:
    type: File
    outputBinding:
      glob: "*.insert.wt.txt"

  spectratype_ndn_wt_file:
    type: File
    outputBinding:
      glob: "*.ndn.wt.txt"

  spectratype_aa_wt_file:
    type: File?
    outputBinding:
      glob: "*.aa.wt.txt"

  spectratype_nt_wt_file:
    type: File?
    outputBinding:
      glob: "*.nt.wt.txt"


baseCommand: ["vdjtools", "CalcSpectratype"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

label: "VDJtools Calc Spectratype"
s:alternateName: "Calculates spectratype, that is, histogram of read counts by CDR3 nucleotide length"

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
  
  Calculates spectratype, that is, histogram of read counts by CDR3 nucleotide length. The
  spectratype is useful to detect pathological and highly clonal repertoires, as the spectratype
  of non-expanded T- and B-cells has a symmetric gaussian-like distribution.
