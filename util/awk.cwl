#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: awk
doc: pattern scanning and processing with awk
$namespaces:
  s: http://schema.org/

requirements:
  DockerRequirement:
    dockerPull: busybox:latest
  InlineJavascriptRequirement: {}

inputs:
  field_separator:
    doc: string that separaters fields on a line
    type: string?
    inputBinding:
      prefix: -F
  inherit_format:
    doc: copy format from target_files input to result output
    type: boolean?
  program:
    doc: AWK program
    type: string
    inputBinding:
      prefix: -e
  target_files:
    doc: files to apply program to
    type: File[]
    inputBinding:
      position: 100
    streamable: true
  variable_setting:
    doc: one or more strings of the format VAR=VALUE to set VAR to VALUE
    type: string[]?
    inputBinding:
      prefix: -v

outputs:
  result:
    type: stdout
    format: |-
      $(inputs.inherit_format && inputs.target_files[0].format ? inputs.target_files[0].format : null)
    streamable: true
stdout: "awk_result$( inputs.inherit_format ? inputs.target_files[0].nameext : '.txt')"

baseCommand: awk

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10
$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf
s:author:
- class: s:Person
  s:email: mailto:pvh@sanbi.ac.za
  s:identifier: https://orcid.org/0000-0001-6553-5274
  s:name: Peter van Heusden
