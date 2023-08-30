#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - var default_output_filename = function() {
            if (inputs.output_filename == ""){
              var root = inputs.segment_file.basename.split('.').slice(0,-1).join('.');
              return (root == "")?inputs.segment_file.basename+".tsv":root+".tsv";
            } else {
              return inputs.output_filename;
            }
          };
hints:
  DockerRequirement:
    dockerPull: biowardrobe2/gat:v0.0.1
  SoftwareRequirement:
    packages:
      gat:
        specs:
          - https://identifiers.org/biotools/gat
          - https://identifiers.org/rrid/RRID:SCR_020949

inputs:

  segment_file:
    type: File
    format: edam:format_3003  # BED
    inputBinding:
      position: 5
      prefix: "-s"
    doc: |
      BED file (strictly 3 columns) with sets of intervals whose association is tested with annotation_file

  annotation_file:
    type: File
    format: edam:format_3003  # BED
    inputBinding:
      position: 6
      prefix: "-a"
    doc: |
      BED file (strictly 3 columns) with sets of intervals that are used for testing association of segment_file

  workspace_file:
    type: File
    format: edam:format_3003  # BED
    inputBinding:
      position: 7
      prefix: "-w"
    doc: |
      BED file (strictly 3 columns) with genomic regions accessible for simulation

  output_filename:
    type: string?
    inputBinding:
      position: 8
      valueFrom: $(default_output_filename())
      prefix: "-S"
    default: ""
    doc: |
      Output report file name

  iterations:
    type: int?
    inputBinding:
      position: 9
      prefix: "-n"
    doc: |
      Number of iterations

  counter:
    type:
    - "null"
    - type: enum
      name: "counter"
      symbols:
      - "nucleotide-overlap"
      - "nucleotide-density"
      - "segment-overlap"
      - "annotation-overlap"
      - "segment-midoverlap"
      - "annotation-midoverlap"
    inputBinding:
      position: 10
      prefix: "-c"
    doc: |
      Set the measure of association that is tested.
      Default: nucleotide-overlap

  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: "-t"
    doc: |
      Threads number

  seed:
    type: int?
    inputBinding:
      position: 12
      prefix: "--random-seed="
      separate: false
    doc: |
      Random seed to initialize random number generator with


outputs:

  report_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Report file

baseCommand: ["gat-run.py", "--ignore-segment-tracks"]

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

s:name: "gat-run"
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
  A common question in genomic analysis is whether two sets of genomic intervals overlap significantly.
  This question arises, for example, in the interpretation of ChIP-Seq or RNA-Seq data. The Genomic
  Association Tester (GAT) is a tool for computing the significance of overlap between multiple sets of
  genomic intervals. GAT estimates significance based on simulation.

  Gat implemements a sampling algorithm. Given a chromosome (workspace) and segments of interest, for
  example from a ChIP-Seq experiment, gat creates randomized version of the segments of interest falling
  into the workspace. These sampled segments are then compared to existing genomic annotations.

  Note:
  --ignore-segment-tracks parameter is hardcoded
