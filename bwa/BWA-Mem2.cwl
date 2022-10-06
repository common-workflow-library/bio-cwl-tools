#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: |
  map medium and long reads (> 100 bp) against reference genome

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_2
  SoftwareRequirement:
    packages:
      bwa-mem2:
        version: [ 2.2.1 ]
        specs: [ https://bio.tools/bwa-mem2 ]

baseCommand: [ bwa-mem2, mem ]

inputs:
  reference_genome:
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles:
     - .bwt.2bit.64
     - .ann
     - .amb
     - .pac
     - ".0123"
    inputBinding:
      position: 1
  reads:
    type: File
    label: "First (forward) set of reads (or unpaired reads)."
    format:
       - edam:format_1929 # FASTA
       - edam:format_1932 # FASTQ-sanger
    inputBinding:
      position: 2
  paired_reads_2:
    type: File?
    label: "Second (reverse) set of reads (if using paired reads)."
    format:
       - edam:format_1929 # FASTA
       - edam:format_1932 # FASTQ-sanger
    inputBinding:
      position: 3
  read_group_header_line:
    type: string?
    doc: |
      read group header line such as '@RG\tID:foo\tSM:bar'
    inputBinding:
      prefix: "-R"

      
arguments:
 - -t
 - $(runtime.cores)
 - -v
 - "1" # Verbosity is set to 1 (errors only)

stdout: |
 ${
    var filename = inputs.reads.nameroot;
    if (inputs.paired_reads_2) {
      filename += "_" + inputs.paired_reads_2.nameroot;
    }
    return filename += ".sam"
  }

outputs:
  aligned_reads:
    type: File
    format: edam:format_2573  # SAM
    outputBinding:
      glob: "*.sam"
    
$namespaces:
  edam: https://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
