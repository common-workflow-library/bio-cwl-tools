#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
id: kraken2
baseCommand:
  - kraken2
inputs:
  database:
    type: 
      - Directory
      - File
    label: "Kraken 2 DB"
    doc: "(either a File refer to the hash.k2d file in the DB or a Directory to reference the entire directory)"
    inputBinding:
      position: 1
      prefix: --db
      valueFrom: |
        ${ return (self.class == "File") ? self.dirname : self.path }
    secondaryFiles:
      - $("opts.k2d")
      - $("taxo.k2d")
  input_sequences:
    type: 
      - File
      - File[]
    label: "Input sequence files"
    format:
      - edam:format_1929  # FASTA
      - edam:format_1930  # FASTQ
    inputBinding:
      position: 2
  output:
    type: string
    label: "Filename for output"
    inputBinding:
      position: 0
      prefix: --output
  threads:
    type: int?
    label: "Number of threads"
    default: 1
    inputBinding:
      position: 0
      prefix: --threads
  quick:
    type: boolean?
    label: "Quick operation (use first hit or hits)"
    inputBinding:
      position: 0
      prefix: --quick
  unclassified_output:
    type: string?
    label: "Print unclassified sequences to this filename"
    inputBinding:
      position: 0
      prefix: unclassified_output
  classified_output:
    type: string?
    label: "Print classified sequences to this filename"
    inputBinding:
      position: 0
      prefix: classified_output
  confidence:
    type: float?
    label: "Confidence score threshold"
    default: 0.0
    inputBinding:
      position: 0
      prefix: --confidence
  minimum-base-quality:
    type: int?
    label: "Minimum base quality used in classification (only used with FASTQ input"
    default: 0
    inputBinding:
      position: 0
      prefix: --minimum-base-quality
  report:
    type:
      - "null"
      - type: record
        name: report_parameters
        fields:
          output_report:
            type: string?
            label: "Print a report with aggregate counts/clade to file"
            inputBinding:
              position: 0
              prefix: --report
          use-mpa-style:
            type: boolean?
            label: "With --report, format report output like Kraken 1's kraken-mpa-report"
            inputBinding:
              position: 0
              prefix: --use-mpa-style
          report-zero-counts:
            type: boolean?
            label: "With --report, report countrs for ALL taxa, even if counts are zero"
            inputBinding:
              position: 0
              prefix: --report-zero-counts
  memory-mapping:
    type: boolean?
    label: "Avoid loading database into RAM"
    inputBinding:
      position: 0
      prefix: --memory-mapping
  paired:
    type: boolean?
    label: "The filenames provided have paired end reads"
    inputBinding:
      position: 0
      prefix: --paired

  use-names:
    type: boolean?
    label: "Print scientific names instead of just taxids"
    inputBinding:
      position: 0
      prefix: --use-names
  gzip-compressed:
    type: boolean?
    label: "Input files are compressed with GZIP"
    inputBinding:
      position: 0
      prefix: --gzip-compressed
  bzip2-compressed:
    type: boolean?
    label: "Input files are compressed with BZIP2"
    inputBinding:
      position: 0
      prefix: --bzip2-compressed
outputs:
  kraken_output:
    type: File
    outputBinding:
      glob: $(inputs.output)
  kraken_report:
    type: File?
    outputBinding:
      glob: $(inputs.output_report)
  classfied_sequences:
    type: File?
    outputBinding:
      glob: $(inputs.classified_output)
  unclassified_sequences:
    type: File?
    outputBinding:
      glob: $(inputs.unclassified_output)
hints:
  - class: SoftwareRequirement
    packages:
      kraken2:
        version:
          - 2.0.8-beta
        specs:
          - http://identifiers.org/biotools/kraken2

requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 45000  # kraken2 standard DB is 38 GB. RAM requirement is modelled on testing with this 
    coresMin: 1
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0"

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
  - "http://edamontology.org/EDAM.owl"
  - "http://schema.org/version/latest/schema.rdf"

s:name: "kraken2"
s:license: "https://spdx.org/licenses/MIT.html"

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "South African National Bioinformatics Institute"
  s:member:
  - class: s:Person
    s:name: Peter van Heusden
    s:email: mailto:pvh@sanbi.ac.za
    s:sameAs:
    - id: https://orcid.org/0000-0001-6553-5274
