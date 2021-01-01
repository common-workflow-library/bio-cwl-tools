#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/

requirements:
  InlineJavascriptRequirement: {}

inputs:
  actg_only:
    doc: only output columns containing exclusively ACTG
    type: boolean?
    inputBinding:
      prefix: -c
  monomorphic_sites:
    doc: output monomorphic sites (suitable for BEAST)
    type: boolean?
    inputBinding:
      prefix: -b
  output_options:
    type:
    - name: output_fasta
      type: record
      fields:
        output_fasta:
          type: boolean?
          inputBinding:
            prefix: -m
    - name: output_phylip
      type: record
      fields:
        output_phylip:
          type: boolean?
          inputBinding:
            prefix: -p
    - name: output_vcf
      type: record
      fields:
        output_vcf:
          type: boolean?
          inputBinding:
            prefix: -v
    default:
      output_fasta: true
  output_pseudo_reference:
    doc: output internal pseudo reference sequence
    type: boolean?
    inputBinding:
      prefix: -r
  sequences:
    doc: FASTA multiple sequence alignment
    type: File
    format: edam:format_1929
    inputBinding:
      position: 100
  site_count_only:
    doc: only output count of constant sites (suitable for IQ-TREE -fconst)
    type: boolean?
    inputBinding:
      prefix: -C

outputs:
  output_sequences:
    type: stdout
    format: |
      ${
        if (inputs.output_options.output_fasta) {
          return 'edam:format_1929';
        } else if (inputs.output_options.output_phylip) {
          return 'edam:format_1997';
        } else if (inputs.output_options.output_vcf) {
          return 'edam:format_3016';
        }
      }
    streamable: true
stdout: |
  ${
    var name = inputs.sequences.nameroot;
    if (inputs.output_options.output_fasta) {
      name += '_snp-sites.fasta';
    } else if (inputs.output_options.output_phylip) {
      name += '_snp-sites.phylip';
    } else if (inputs.output_options.output_vcf) {
      name += '_snp-sites.vcf';
    }
    return name;
  }

baseCommand: snp-sites

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0
  SoftwareRequirement:
    packages:
      snp-sites:
        specs:
        - https://anaconda.org/bioconda/snp-sites
      versions:
      - 2.5.1
$schemas:
- https://schema.org/version/latest/schemaorg-current-http.rdf
- http://edamontology.org/EDAM_1.18.owl
s:author:
- class: s:Person
  s:email: mailto:pvh@sanbi.ac.za
  s:identifier: https://orcid.org/0000-0001-6553-5274
  s:name: Peter van Heusden
