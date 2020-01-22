#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2

hints:
  SoftwareRequirement:
    packages:
      bismark:
        specs: [ "http://identifiers.org/biotools/bismark" ]
        version: [ "0.0.2" ]

inputs:

  alignment_report:
    type: File
    inputBinding:
      prefix: "--alignment_report"
    label: "Bismark alignment and methylation report"
    doc: "Bismark generated alignment and methylation summary report"

  splitting_report:
    type: File?
    inputBinding:
      prefix: "--splitting_report"
    label: "Methylation extraction log"
    doc: "Log file giving summary statistics about methylation extraction"

  mbias_report:
    type: File?
    inputBinding:
      prefix: "--mbias_report"
    label: "Methylation bias plot"
    doc: "QC data showing methylation bias across read lengths"

  deduplication_report:
    type: File?
    inputBinding:
      prefix: "--dedup_report"
    label: "Bismark deduplication report"
    doc: "Bismark generated deduplication report"

  nucleotide_report:
    type: File?
    inputBinding:
      prefix: "--nucleotide_report"
    label: "Nucleotide coverage report"
    doc: "Bismark nucleotide coverage report"

outputs:

  collected_report:
    type: File
    label: "HTML report page"
    doc: "Bismark generated graphical HTML report page"
    outputBinding:
      glob: "*"

baseCommand: ["bismark2report"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bismark_report"
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
  This tool uses a Bismark alignment report to generate a graphical HTML report page. Optionally, further reports of
  the Bismark suite such as deduplication, methylation extractor splitting or M-bias reports can be specified as well.

  Skipped arguments:
    -o/--output
    --dir