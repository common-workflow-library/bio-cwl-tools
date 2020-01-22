#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie:v1.2.0
- class: SoftwareRequirement
  packages:
    bowtie:
      specs: [ "http://identifiers.org/biotools/bowtie" ]
      version: [ "1.2.0" ]

inputs:

  fasta_file:
    type:
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 25
    doc: |
      comma-separated list of files with ref sequences

  index_base_name:
    type: string
    inputBinding:
      position: 26
    doc: |
      write Ebwt data to files with this dir/basename

  force_large_index:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 3
      prefix: '--large-index'
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides

  color:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 4
      prefix: '--color'
    doc: |
      build a colorspace index

  noauto:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '--noauto'
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting

  packed:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 6
      prefix: '--packed'
    doc: |
      use packed strings internally; slower, less memory

  bmax:
    type:
      - "null"
      - int
    inputBinding:
      position: 7
      prefix: '--bmax'
    doc: |
      max bucket sz for blockwise suffix-array builder

  bmaxdivn:
    type:
      - "null"
      - int
    inputBinding:
      position: 8
      prefix: '--bmaxdivn'
    doc: |
      max bucket sz as divisor of ref len (default: 4)

  dcv:
    type:
      - "null"
      - int
    inputBinding:
      position: 9
      prefix: '--dcv'
    doc: |
      diff-cover period for blockwise (default: 1024)

  nodc:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 10
      prefix: '--nodc'
    doc: |
      disable diff-cover (algorithm becomes quadratic)

  noref:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 11
      prefix: '--noref'
    doc: |
      don't build .3/.4 index files

  justref:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 12
      prefix: '--justref'
    doc: |
      just build .3/.4 index files

  offrate:
    type:
      - "null"
      - int
    inputBinding:
      position: 13
      prefix: '--offrate'
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)

  ftabchars:
    type:
      - "null"
      - int
    inputBinding:
      position: 14
      prefix: '--ftabchars'
    doc: |
      # of chars consumed in initial lookup (default: 10)

  ntoa:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 15
      prefix: '--ntoa'
    doc: |
      convert Ns in reference to As

  seed:
    type:
      - "null"
      - int
    inputBinding:
      position: 16
      prefix: '--seed'
    doc: |
      seed for random number generator

  quiet:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 17
      prefix: '--quiet'
    doc: |
      verbose output (for debugging)

outputs:

  indices:
    type: File[]
    outputBinding:
      glob: "*"
      outputEval: |
        ${
          var output_array = [];
          for (var i = 0; i < self.length; i++){
            if (self[i].class == "File"){
              output_array.push(self[i]);
            }
          }
          return output_array;
        }

baseCommand:
  - bowtie-build

arguments:
  - valueFrom: $('> ' + inputs.index_base_name + '.log')
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bowtie_build"
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
  Tool runs bowtie-build
  Not supported parameters:
    -c  -  reference sequences given on cmd line (as <seq_in>)