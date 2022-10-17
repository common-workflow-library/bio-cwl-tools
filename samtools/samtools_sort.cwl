#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

doc: Sort a bam file by read names.

requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 15000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.14--hb421002_0

baseCommand: ["samtools", "sort"]
arguments:
  - valueFrom: $(runtime.cores)
    prefix: --threads  # a.k.a -@
  - prefix: -m
    valueFrom: ${ return(parseInt(runtime.ram/runtime.cores-100).toString() + "M") }
    position: 1
    # specifies the allowed maximal memory usage per thread before
    # samtools start to outsource memory to temporary files
  - prefix: -T
    valueFrom: $(runtime.tmpdir)
  - --no-PG  # don't put the filepaths in the metadata, makes the output more reproducible

inputs:
  unsorted_alignments:
    doc: aligned reads to be checked in sam or bam format
    type: File
    format:
      - edam:format_2572  # BAM
      - edam:format_2573  # SAM
    inputBinding:
      position: 2
  by_name:
    doc: If true, will sort by name, otherwise will sort by genomic position
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -n
  force_format:
    doc: If true, will force binary output (BAM)
    type:
      - 'null'
      - type: enum
        symbols: [ SAM, BAM, CRAM ]
    inputBinding:
      prefix: -O

stdout: |
  ${
    var filename = inputs.unsorted_alignments.nameroot + "_sorted.";
    if (inputs.force_format !== undefined) {
      return filename + inputs.force_format.toLowerCase();
    }
    return filename + (inputs.unsorted_alignments.format == "http://edamontology.org/format_2572" ? "bam" : "sam");
   }

outputs:
  sorted_alignments:
    type: stdout
    format: |
      ${
         if (inputs.force_format === undefined) {
           return inputs.unsorted_alignments.format;
         }
         if (inputs.force_format == "SAM") {
           return "http://edamontology.org/format_2573";
         }
         if (inputs.force_format == "BAM") {
           return "http://edamontology.org/format_2572";
         }
         if (inputs.force_format == "CRAM") {
           return "http://edamontology.org/format_3462";
         }
       }

$namespaces:
  edam: http://edamontology.org/ 
