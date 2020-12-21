#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: Assign Nextstrain clades to SARS-CoV-2 sequences and provide QC information
id: nextclade
label: Nextclade

s:author:
  - class: s:Organization
    s:name: Neher Lab
    s:url: https://neherlab.org/

s:contributor:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-6553-5274
    s:email: mailto:pvh@sanbi.ac.za
    s:name: Peter van Heusden

s:codeRepository: https://github.com/nextstrain/nextclade

requirements:
  DockerRequirement:
    dockerPull: neherlab/nextclade:0.10.0-alpine

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 512  # 512 MB

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
    format: iana:application/json
    outputBinding:
      glob: $(inputs.output_options.output_json_filename)
  output_csv:
    type: File?
    format: iana:text/csv  # Comma-separated values
    outputBinding:
      glob: $(inputs.output_options.output_csv_filename)
  output_tsv_clades_only:
    type: File?
    format: iana:text/tab-separated-values  # Tab-separated values
    outputBinding:
      glob: $(inputs.output_options.output_tsv_clades_only_filename)
  output_tsv:
    type: File?
    format: iana:text/tab-separated-values
    outputBinding:
      glob: $(inputs.output_options.output_tsv_filename)
  output_tree:
    type: File?
    format: iana:application/json
    outputBinding:
      glob: $(inputs.output_options.output_tree_filename)

baseCommand: [ nextclade.js ]

  
$namespaces:
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf
  - http://edamontology.org/EDAM_1.18.owl

