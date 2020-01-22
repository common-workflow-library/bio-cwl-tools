#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_output_filename = function(ext) {
        var alt_ext = "";
        if (inputs.input_file_type == "bam") {
          alt_ext = ".sorted.bam";
        } else if (inputs.input_file_type == "bigwig") {
          alt_ext = ".bw";
        } else if (inputs.input_file_type == "bed") {
          ext = ".bedGraph";
        } else {
          alt_ext = "";
        }
        ext = (ext || ext=="")?ext:alt_ext;
        if (inputs.output_basename == ""){
          var root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.input_file.basename+ext:root+ext;
        } else {
          return inputs.output_basename+ext;
        }
      };
    - var get_log_filename = function() {
        var ext = ".log";
        if (inputs.output_basename == ""){
          var root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.input_file.basename+ext:root+ext;
        } else {
          return inputs.output_basename+ext;
        }
      };

hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/crossmap:0.2.7--py27_0

inputs:

  input_file_type:
    type: string
    inputBinding:
      position: 1
    doc: |
      bam	    convert alignment file in BAM format.
      bed	    convert genome cooridnate or annotation file in BED or BED-like format.
      bigwig	convert genome coordinate file in BigWig format.

  chain_file:
    type: File
    inputBinding:
      position: 2
    doc: |
      Chain file

  input_file:
    type: File
    inputBinding:
      position: 3
    secondaryFiles: |
      ${
        return (inputs.input_file_type == "bam")?self.basename+".bai":[];
      }
    doc: |
      Input file BAM(+bai), BED, BigWig.

  output_basename:
    type: string?
    inputBinding:
      position: 4
      valueFrom: $(get_output_filename(""))
    default: ""
    doc: |
      Name for the generated output file

  bam_insert_size:
    type: int?
    inputBinding:
      position: 5
      prefix: -m
    doc: |
      For BAM only: Average insert size of pair-end sequencing (bp). [default=200.0]

  bam_stdev:
    type: float?
    inputBinding:
      position: 6
      prefix: -s
    doc: |
      For BAM only: Stanadard deviation of insert size. [default=30.0]

  bam_fold:
    type: float?
    inputBinding:
      position: 7
      prefix: -t
    doc: |
      For BAM only: A mapped pair is considered as "proper pair" if both
                    ends mapped to different strand and the distance
                    between them is less then '-t' * stdev from the mean.
                    [default=3.0]

  bam_append_tags:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -a
    doc: |
      For BAM only: Add tag to each alignment

outputs:

  projected_file:
    type: File
    outputBinding:
      glob: |
        ${
          return get_output_filename();
        }
    secondaryFiles: |
      ${
        if (inputs.input_file_type == "bam") {
          return self.basename + ".bai";
        } else {
          return "null";
        }
      }
    doc: |
      Projected output file

  unmap_file:
    type: File?
    outputBinding:
      glob: "*.unmap"
    doc: |
      Unmap output file

  log_file:
    type: File
    outputBinding:
      glob: $(get_log_filename())
    doc: |
      Log file

stdout: $(get_log_filename())

baseCommand: ["CrossMap.py"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "crossmap"
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
  Runs CrossMap.py script to project input BAM, BED, BIGWIG file based on input chain file.
  Not supported input file types: SAM, GFF, VCF, WIG

  If `output_basename` is not set, call get_output_filename() and get_log_filename() functions to
  get default output and log filenames. Input `output_basename` should not include extension.