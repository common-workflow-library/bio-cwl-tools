#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Assign Nextstrain clades to SARS-CoV-2 sequences and provide QC information
id: nextclade
label: Nextclade

dct:creator:
  "@id": "https://orcid.org/0000-0001-6553-5274"
  foaf:name: Peter van Heusden
  foaf:mbox: "mailto:pvh@sanbi.ac.za"

requirements:
  DockerRequirement:
    dockerPull: neherlab/nextclade:0.8.1-alpine

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 512  # 512 MB

inputs:
  input_fasta:
    type: File
    doc: .fasta or .txt file with input sequences
    format: 
      - edam:format_1929  # FASTA
      - edam:format_1964  # plain text format
    inputBinding:
      prefix: --input-fasta
  input_qc_config:
    type: File?
    doc: QC config json file containing custom QC configuration
    format: edam:format_3464  # JSON
    inputBinding:
      prefix: --input-qc-config
  input_root_seq:
    type: File?
    doc: plain text file containing custom root sequence
    format: edam:format_1964
    inputBinding:
      prefix: --input-root-seq
  input-tree:
    type: File?
    doc: Auspice JSON v2 file containing custom reference tree
    format: edam:format_3464
    inputBinding:
      prefix: --input-tree
  input-gene-map:
    type: File?
    doc: 'JSON file containing custom gene map. Gene map (sometimes also called "gene annotations") is used to resolve aminoacid changes in genes.'
    format: edam:format_3464
    inputBinding:
      prefix: --input-gene-map
  input-pcr-primers:
    type: File?
    doc: CSV file containing a list of custom PCR primer sites. These are used to report mutations in these sites.
    format: edam:format_3572
    inputBinding:
      prefix: --input-pcr-primers
  output_options:
    type:
      type: record
      name: output_options_record
      fields:
        output_json_filename:
          type: string?
          doc: Filename of output JSON results file
          inputBinding:
            prefix: --output-json
        output_csv_filename:
          type: string?
          doc: Filename of output CSV results file
          inputBinding:
            prefix: --output-csv
        output_tsv_filename:
          type: string?
          doc: Filename of output TSV results file
          inputBinding:
            prefix: --output-tsv
        output_tsv_clades_only_filename:
          type: string?
          doc: Filename to output CSV clades-only file
          inputBinding:
            prefix: --output-tsv-clades-only
        output_tree_filename:
          type: string?
          doc: Filename of output Auspice v2 tree file
          inputBinding:
            prefix: --output-tree

outputs:
  output_json:
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: $(inputs.output_options.output_json_filename)
  output_csv:
    type: File?
    format: edam:format_3572  # Comma-separated values
    outputBinding:
      glob: $(inputs.output_options.output_csv_filename)
  output_tsv_clades_only:
    type: File?
    format: edam:format_3475  # Tab-separated values
    outputBinding:
      glob: $(inputs.output_options.output_tsv_clades_only_filename)
  output_tsv:
    type: File?
    format: edam:format_3475
    outputBinding:
      glob: $(inputs.output_options.output_tsv_filename)
  output_tree:
    type: File?
    format: edam:format_3464
    outputBinding:
      glob: $(inputs.output_options.output_tree_filename)

baseCommand: [ nextclade.js ]

  
$namespaces:
  edam: http://edamontology.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl

