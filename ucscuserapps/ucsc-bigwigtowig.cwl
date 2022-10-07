#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - var default_output_filename = function() {
            var basename = inputs.bigwig_file.location.split('/').slice(-1)[0];
            var root = basename.split('.').slice(0,-1).join('.');
            var ext = ".wig";
            return (root == "")?basename+ext:root+ext;
          };
hints:
  DockerRequirement:
    dockerPull: biowardrobe2/ucscuserapps:v358

inputs:

  bigwig_file:
    type: File
    inputBinding:
      position: 1
    doc: |
      Input bigWig file

  chrom:
    type: string?
    inputBinding:
      position: 2
      prefix: -chrom=
      separate: false
    doc: |
      if set restrict output to given chromosome

  start_pos:
    type: int?
    inputBinding:
      position: 3
      prefix: -start=
      separate: false
    doc: |
      if set, restrict output to only that over start

  end_pos:
    type: int?
    inputBinding:
      position: 4
      prefix: -end=
      separate: false
    doc: |
      if set, restict output to only that under end

  output_filename:
    type: string?
    inputBinding:
      position: 5
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      If set, writes the output wig file to output_filename,
      otherwise generates filename from default_output_filename()


outputs:
  wig_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == ""){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }

baseCommand: bigWigToWig

$namespaces:
  s: http://schema.org/
  edam: https://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "ucsc-bigwigtowig"
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
  Tool converts bigWig to Wig file. If input bigWig file was generated from bedGraph, the tool will
  return output in bedGraph format.

  `default_output_filename` function returns filename for generated Wig if `output_filename` is not provided.
  Default filename is generated on the base of `bigwig_file` basename with the updated to `*.wig` extension.
