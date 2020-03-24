#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sratoolkit:v2.8.2-1
- class:  SoftwareRequirement
  packages:
    sra-tools:
      specs: [ "http://identifiers.org/biotools/sra-tools" ]
      version: [ "2.8.2" ]

inputs:

  sra_file:
    type: File
    inputBinding:
      position: 60
    doc: |
      Input file

  split_spot:
    type: boolean?
    inputBinding:
      position: 2
      prefix: "--split-spot"
    doc: |
      Split spots into individual reads

  min_spot_id:
    type: string?
    inputBinding:
      position: 3
      prefix: "--minSpotId"
    doc: |
      Minimum spot id

  max_spot_id:
    type: string?
    inputBinding:
      position: 4
      prefix: "--maxSpotId"
    doc: |
      Maximum spot id

  clip:
    type: boolean?
    inputBinding:
      position: 6
      prefix: "--clip"
    doc: |
      Clip adapter sequences

  min_read_len:
    type: int?
    inputBinding:
      position: 7
      prefix: "--minReadLen"
    doc: |
      Filter by sequence length >= <len>

  qual_filter:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "--qual-filter"
    doc: |
      Filter used in early 1000 Genomes data: no sequences starting or ending with >= 10N

  qual_filter_1:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--qual-filter-1"
    doc: |
      Filter used in current 1000 Genomes data

  aligned:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--aligned"
    doc: |
      Dump only aligned sequences

  unaligned:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "--unaligned"
    doc: |
      Dump only unaligned sequences

  aligned_region:
    type: string?
    inputBinding:
      position: 13
      prefix: "--aligned-region"
    doc: |
      Filter by position on genome.
      Name can either be accession.version
      (ex:NC_000001.10) or file specific name
      (ex:"chr1" or "1"). "from" and "to" are 1-based coordinates

  matepair_distance:
    type:
      - "null"
      - type: enum
        name: "distance"
        symbols: ["from-to","unknown"]
    inputBinding:
      position: 14
      prefix: "--matepair-distance"
    doc: |
      Filter by distance beiween matepairs.
      Use "unknown" to find matepairs split
      between the references. Use from-to to limit
      matepair distance on the same reference

  skip_technical:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "--skip-technical"
    doc: |
      Dump only biological reads

  split_files:
    type: boolean?
    inputBinding:
      position: 20
      prefix: "--split-files"
    doc: |
      Dump each read into separate file.
      Files will receive suffix corresponding to read number

  split_3:
    type: boolean?
    inputBinding:
      position: 21
      prefix: "--split-3"
    doc: |
      Legacy 3-file splitting for mate-pairs:
      First biological reads satisfying dumping
      conditions are placed in files *_1.fastq and
      *_2.fastq If only one biological read is
      present it is placed in *.fastq Biological
      reads and above are ignored.

  dumpcs:
    type: string?
    inputBinding:
      position: 25
      prefix: "--dumpcs"
    doc: |
      Formats sequence using color space (default
      for SOLiD),"cskey" may be specified for
      translation

  dumpbase:
    type: boolean?
    inputBinding:
      position: 26
      prefix: "--dumpbase"
    doc: |
      Formats sequence using base space (default
      for other than SOLiD).

  offset:
    type: int?
    inputBinding:
      position: 27
      prefix: "--offset"
    doc: |
      Offset to use for quality conversion, default is 33

  fasta:
    type: int?
    inputBinding:
      position: 28
      prefix: "--fasta"
    doc: |
      FASTA only, no qualities, optional line
      wrap width (set to zero for no wrapping)

  suppress_qual_for_cskey:
    type: boolean?
    inputBinding:
      position: 29
      prefix: "--suppress-qual-for-cskey"
    doc: |
      supress quality-value for cskey

  origfmt:
    type: boolean?
    inputBinding:
      position: 30
      prefix: "--origfmt"
    doc: |
      Defline contains only original sequence name

  readids:
    type: boolean?
    inputBinding:
      position: 31
      prefix: "--readids"
    doc: |
      Append read id after spot id as 'accession.spot.readid' on defline

  helicos:
    type: boolean?
    inputBinding:
      position: 32
      prefix: "--helicos"
    doc: |
      Helicos style defline

  defline_seq:
    type: string?
    inputBinding:
      position: 33
      prefix: "--defline-seq"
    doc: |
      Defline format specification for sequence.

  defline_qual:
    type: string?
    inputBinding:
      position: 34
      prefix: "--defline-qual"
    doc: |
      Defline format specification for quality.
      <fmt> is string of characters and/or
      variables. The variables can be one of: $ac
      - accession, $si spot id, $sn spot
      name, $sg spot group (barcode), $sl spot
      length in bases, $ri read number, $rn
      read name, $rl read length in bases. '[]'
      could be used for an optional output: if
      all vars in [] yield empty values whole
      group is not printed. Empty value is empty
      string or for numeric variables. Ex:
      @$sn[_$rn]/$ri '_$rn' is omitted if name
      is empty

  disable_multithreading:
    type: boolean?
    inputBinding:
      position: 35
      prefix: "--disable-multithreading"
    doc: |
      disable multithreading

outputs:
  fastq_file_1:
    type: File
    outputBinding:
      glob: |
        ${
          return [
              inputs.sra_file.basename.split(".")[0] + ".fastq",
              inputs.sra_file.basename.split(".")[0] + "_1.fastq"
            ];
        }

  fastq_file_2:
    type: File?
    outputBinding:
      glob: |
        ${
          return inputs.sra_file.basename.split(".")[0] + "_2.fastq";
        }
  all_fastq_files:
    type: File[]
    outputBinding:
      glob: |
        ${
          return inputs.sra_file.basename.split(".")[0] + "*.fastq";
        }

baseCommand: [fastq-dump]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "fastq_dump"
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
  Tool runs fastq-dump from NCBI SRA toolkit
  Supports only file inputs.
  Output file names are formed on the base of `sra_file` input basename.