#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.bed_file.basename;
          } else {
            return inputs.output_filename;
          }
        };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
- class: SoftwareRequirement
  packages:
    bedtools:
      specs: [ "http://identifiers.org/biotools/bedtools" ]
      version: [ "2.26.0" ]

inputs:

  bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "The input BED file must be sorted by chrom, then start"

  max_distance:
    type: int?
    inputBinding:
      position: 6
      prefix: "-d"
    doc: "Maximum distance between features to be merged" 

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"

outputs:

  merged_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Merged BED file"

baseCommand: ["bedtools", "merge"]
stdout: $(default_output_filename())

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bedtools_merge"
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
  Merges features from BED file. Only selected parameters are implemented.