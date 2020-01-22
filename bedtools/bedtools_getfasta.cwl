#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.intervals_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.intervals_file.basename+".fa":root+".fa";
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

  genome_fasta_file:
    type: File
    secondaryFiles: $(self.basename+".fai")  # due to bug in cwltool==1.0.20190621234233
    inputBinding:
      position: 5
      prefix: "-fi"
    doc: "Genome file in FASTA format"

  intervals_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-bed"
    doc: "Intervals file defined in a BED/GFF/VCF format"

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"

outputs:

  sequences_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Sequences file"

baseCommand: ["bedtools", "getfasta"]
stdout: $(default_output_filename())

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bedtools_getfasta"
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
  Extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file. Only selected parameters are implemented.