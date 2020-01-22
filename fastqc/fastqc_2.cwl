#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqc:v0.11.5
- class:  SoftwareRequirement
  packages:
    fastqc:
      specs: [ "http://identifiers.org/biotools/fastqc" ]
      version: [ "0.1.11.5" ]

inputs:

  reads_file:
    type:
      - File
    inputBinding:
      position: 50
    doc: |
      Input bam,sam,bam_mapped,sam_mapped or fastq file

  format_enum:
    type:
      - "null"
      - type: enum
        name: "format"
        symbols: ['bam','sam','bam_mapped','sam_mapped','fastq']
    inputBinding:
      position: 6
      prefix: '--format'
    doc: |
      Bypasses the normal sequence file format detection and
      forces the program to use the specified format.  Valid
      formats are bam,sam,bam_mapped,sam_mapped and fastq

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: '--threads'
    doc: |
      Specifies the number of files which can be processed
      simultaneously.  Each thread will be allocated 250MB of
      memory so you shouldn't run more threads than your
      available memory will cope with, and not more than
      6 threads on a 32 bit machine

  contaminants:
    type:
      - "null"
      - File
    inputBinding:
      position: 8
      prefix: '--contaminants'
    doc: |
      Specifies a non-default file which contains the list of
      contaminants to screen overrepresented sequences against.
      The file must contain sets of named contaminants in the
      form name[tab]sequence.  Lines prefixed with a hash will
      be ignored.

  adapters:
    type:
      - "null"
      - File
    inputBinding:
      position: 9
      prefix: '--adapters'
    doc: |
      Specifies a non-default file which contains the list of
      adapter sequences which will be explicity searched against
      the library. The file must contain sets of named adapters
      in the form name[tab]sequence.  Lines prefixed with a hash
      will be ignored.

  limits:
    type:
      - "null"
      - File
    inputBinding:
      position: 10
      prefix: '--limits'
    doc: |
      Specifies a non-default file which contains a set of criteria
      which will be used to determine the warn/error limits for the
      various modules.  This file can also be used to selectively
      remove some modules from the output all together.  The format
      needs to mirror the default limits.txt file found in the
      Configuration folder.

  kmers:
    type:
      - "null"
      - int
    inputBinding:
      position: 11
      prefix: '--kmers'
    doc: |
      Specifies the length of Kmer to look for in the Kmer content
      module. Specified Kmer length must be between 2 and 10. Default
      length is 7 if not specified.

  casava:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 13
      prefix: '--casava'
    doc: |
      Files come from raw casava output. Files in the same sample
      group (differing only by the group number) will be analysed
      as a set rather than individually. Sequences with the filter
      flag set in the header will be excluded from the analysis.
      Files must have the same names given to them by casava
      (including being gzipped and ending with .gz) otherwise they
      won't be grouped together correctly.

  nofilter:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 14
      prefix: '--nofilter'
    doc: |
      If running with --casava then don't remove read flagged by
      casava as poor quality when performing the QC analysis.

  hide_group:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 15
      prefix: '--nogroup'
    doc: |
      Disable grouping of bases for reads >50bp. All reports will
      show data for every base in the read.  WARNING: Using this
      option will cause fastqc to crash and burn if you use it on
      really long reads, and your plots may end up a ridiculous size.
      You have been warned!

outputs:

  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
  summary_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }

baseCommand: [fastqc, --extract, --outdir, .]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "fastqc_2"
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool runs FastQC from Babraham Bioinformatics