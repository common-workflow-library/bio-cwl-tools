#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: awk
doc: pattern scanning and processing with awk
$namespaces:
  s: http://schema.org/

requirements:
  DockerRequirement:
    dockerPull: busybox:latest

inputs:
  field_separator:
    doc: string that separaters fields on a line
    type: string?
    inputBinding:
      prefix: -F
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
    streamable: true
stdout: awk_result$(inputs.target_files[0].nameext)

baseCommand: awk

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10
$schemas:
- https://schema.org/version/latest/schemaorg-current-http.rdf
s:author:
- class: s:Person
  s:email: mailto:pvh@sanbi.ac.za
  s:identifier: https://orcid.org/0000-0001-6553-5274
  s:name: Peter van Heusden
