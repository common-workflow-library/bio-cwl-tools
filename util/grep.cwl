#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: grep
doc: print lines that match patterns
$namespaces:
  s: http://schema.org/

requirements:
  DockerRequirement:
    dockerPull: busybox:latest
  InlineJavascriptRequirement: {}

inputs:
  inherit_format:
    doc: copy format from search_target input to result output
    type: boolean?
  ignore_case:
    doc: ignore case while matching
    type: boolean?
    inputBinding:
      prefix: -i
  invert_selection:
    doc: print non-matching lines
    type: boolean?
    inputBinding:
      prefix: -v
  line_no:
    doc: add line number prefix to results
    type: boolean?
    inputBinding:
      prefix: -n
  lines_after:
    doc: print up to N lines of trailing content
    type: int?
    inputBinding:
      prefix: -A
  lines_before:
    doc: print up to N lines of leading content
    type: int?
    inputBinding:
      prefix: -B
  literal_pattern:
    doc: interpret the pattern as a literal string, not a regular expression
    type: boolean?
    inputBinding:
      prefix: -F
  pattern:
    doc: pattern to search for
    type: string
    inputBinding:
      prefix: -e
  search_target:
    doc: file to search
    type: File
    inputBinding:
      position: 100
    streamable: true

outputs:
  result:
    type: stdout
    streamable: true
    format: $(inputs.inherit_format ? inputs.search_target.format : null)

stdout: search_result$(inputs.search_target.nameext)

baseCommand: grep

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
