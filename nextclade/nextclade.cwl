#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Assign Nextstrain clades to SARS-CoV-2 sequences and provide QC information
label: Nextclade

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-6553-5274
    s:email: mailto:pvh@sanbi.ac.za
    s:name: Peter van Heusden

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/nextclade_js:0.14.3--h9ee0642_0
  ResourceRequirement:
    coresMin: 1
    ramMin: 512  # 512 MB
  SoftwareRequirement:
    packages:
      nextclade: 
       version: 
        - 0.14.3
       specs: 
        - https://anaconda.org/bioconda/nextclade_js
        - https://github.com/nextstrain/nextclade

inputs:
  sequences:
    type: File
    doc: .fasta or .txt file with input sequences
    format: 
      - edam:format_1929  # FASTA
      - iana:text/plain  # plain text format
    inputBinding:
      prefix: --input-fasta
  qc_config:
    type: File?
    doc: QC config json file containing custom QC configuration
    format: iana:application/json  # JSON
    inputBinding:
      prefix: --input-qc-config
  root_seq:
    type: File?
    doc: plain text file containing custom root sequence
    format: iana:text/plain
    inputBinding:
      prefix: --input-root-seq
  tree:
    type: File?
    doc: Auspice JSON v2 file containing custom reference tree
    format: iana:application/json
    inputBinding:
      prefix: --input-tree
  gene_map:
    type: File?
    doc: 'JSON file containing custom gene map. Gene map (sometimes also called "gene annotations") is used to resolve aminoacid changes in genes.'
    format: iana:application/json
    inputBinding:
      prefix: --input-gene-map
  pcr_primers:
    type: File?
    doc: CSV file containing a list of custom PCR primer sites. These are used to report mutations in these sites.
    format: iana:text/csv
    inputBinding:
      prefix: --input-pcr-primers
  output_options:
    type:
      type: record
      name: output_options_record
      fields:
        report_json:
          type: boolean?
          doc: Filename of output JSON results file
          inputBinding:
            prefix: --output-json
            valueFrom: $(inputs.sequences.nameroot)_report.json
        report_csv:
          type: boolean?
          doc: Filename of output CSV results file
          inputBinding:
            prefix: --output-csv
            valueFrom: $(inputs.sequences.nameroot)_report.csv
        report_tsv:
          type: boolean?
          doc: Filename of output TSV results file
          inputBinding:
            prefix: --output-tsv
            valueFrom: $(inputs.sequences.nameroot)_report.tsv
        report_tsv_clades_only_filename:
          type: boolean?
          doc: Filename to output CSV clades-only file
          inputBinding:
            prefix: --output-tsv-clades-only
            valueFrom: $(inputs.sequences.nameroot)_clades_report.csv
        tree_filename:
          type: boolean?
          doc: Filename of output Auspice v2 tree file
          inputBinding:
            prefix: --output-tree
            valueFrom: $(inputs.sequences.nameroot)_tree.json
    default:
      report_tsv: true
      report_json: true

outputs:
  report_json:
    type: File?
    format: iana:application/json
    outputBinding:
      glob: $(inputs.sequences.nameroot)_report.json
  report_csv:
    type: File?
    format: iana:text/csv  # Comma-separated values
    outputBinding:
      glob: $(inputs.sequences.nameroot)_report.csv
  report_tsv_clades_only:
    type: File?
    format: iana:text/tab-separated-values  # Tab-separated values
    outputBinding:
      glob: $(inputs.sequences.nameroot)_clades_report.tsv
  report_tsv:
    type: File?
    format: iana:text/tab-separated-values
    outputBinding:
      glob: $(inputs.sequences.nameroot)_report.tsv
  report_tree:
    type: File?
    format: iana:application/json
    outputBinding:
      glob: $(inputs.sequences.nameroot)_tree.json

baseCommand: [ nextclade.js ]

  
$namespaces:
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf
  - https://edamontology.org/EDAM_1.18.owl

