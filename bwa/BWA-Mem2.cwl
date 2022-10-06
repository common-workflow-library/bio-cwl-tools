#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow

label: |
  map medium and long reads (> 100 bp) against reference genome

inputs:
  reference_genome:
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles:
     - .bwt.2bit.64?
     - .ann?
     - .amb?
     - .pac?
     - ".0123?"
  paired_reads_1:
    type: File
    label: "First (forward) set of reads"
    format:
       - edam:format_1929 # FASTA
       - edam:format_1932 # FASTQ-sanger
  paired_reads_2:
    type: File
    label: "Second (reverse) set of reads"
    format:
       - edam:format_1929 # FASTA
       - edam:format_1932 # FASTQ-sanger
  do_auto_name:
    type: boolean
    default: False
    label: "Auto-assign read groups"
    doc: "If true, use the file name to automatically assign the read groups value."

outputs:
  sorted_alignments:
    type: File
    format: edam:format_2572  # BAM
    outputSource: sort/sorted_alignments

steps:
  index_genome:
    run: BWA-Mem2-index.cwl
    in:
      sequences: reference_genome
    when: |
      $(inputs.sequences.secondaryFiles !== undefined)
    out: [ indexed_sequences ]
  compute_read_group_header:
    run: ReadGroup.cwl
    when: $(inputs.do_auto_name)
    in:
      do_auto_name: do_auto_name
      input1: paired_reads_1
      input2: paired_reads_2
    out: [ read_group_name ]
  align:
    run: BWA-Mem2-paired.cwl
    in:
      reference_genome:
        source: [ index_genome/indexed_sequences, reference_genome ]
        pickValue: first_non_null
      paired_reads_1: paired_reads_1
      paired_reads_2: paired_reads_2
      read_group_header_line: compute_read_group_header/read_group_name
    out: [ aligned_reads ]
  sort:
    run: ../samtools/samtools_sort.cwl
    in: 
      unsorted_alignments: align/aligned_reads
      force_format:
        default: BAM
    out: [ sorted_alignments ]

requirements:
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
 
$namespaces:
  edam: https://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
