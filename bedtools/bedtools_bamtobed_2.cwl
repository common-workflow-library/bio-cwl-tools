#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.bam_file.basename+".bed":root+".bed";
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

  bam_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "Input BAM file (not necessary sorted or indexed)"

  output_filename:
    type: string?
    default: ""
    doc: "Output BED filename"

outputs:

  bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Sequences file"

baseCommand: ["bedtools", "bamtobed"]
stdout: $(default_output_filename())

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bedtools_bamtobed_2"
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
  Converts BAM to BED. All additional options are not implemented.