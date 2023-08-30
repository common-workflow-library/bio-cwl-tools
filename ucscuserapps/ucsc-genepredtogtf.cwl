#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358_2


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.annotation_tsv_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.annotation_tsv_file.basename+".gtf":root+".gtf";
          } else {
            return inputs.output_filename;
          }
        };


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      TSV_FILE=$0
      GTF_FILE=$1
      cut -f 2- $TSV_FILE | grep -v "exonCount" | genePredToGtf file stdin $GTF_FILE
    inputBinding:
      position: 5
    doc: |
      Bash function to run cut -f 2- in.gp | genePredToGtf file stdin out.gp

  annotation_tsv_file:
    type: File
    inputBinding:
      position: 6
    doc: "TSV annotation file from UCSC Genome Browser"

  output_filename:
    type: string?
    default: ""
    inputBinding:
      valueFrom: $(default_output_filename())
      position: 7
    doc: "Output file name"


outputs:

  annotation_gtf_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "GTF annotation file"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]


stdout: genepredtogtf_stdout.log
stderr: genepredtogtf_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

s:name: "ucsc-genepredtogtf"
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
  genePredToGtf - Convert genePred table or file to gtf
