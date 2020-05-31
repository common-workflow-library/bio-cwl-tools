#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie2:v2.3.0
  dockerFile: >
    $import: ./dockerfiles/bowtie2-Dockerfile
- class: SoftwareRequirement
  packages:
      bowtie2:
        specs: [ "http://identifiers.org/biotools/bowtie2" ]
        version: [ "2.3.0" ]

inputs:

  reference_in:
    type:
      - File
      - type: array
        items: File
    doc: >
      comma-separated list of files with ref sequences
    inputBinding:
      itemSeparator: ","
      position: 25

  bt2_index_base:
    type: string?
    doc: |
      write bt2 data to files with this dir/basename
    inputBinding:
      position: 26
      shellQuote: false
    default: bowtie2_indices

  f:
    type:
      - "null"
      - boolean
    doc: |
      reference files are Fasta (default)
    inputBinding:
      position: 1
      prefix: '-f'

  c:
    type:
      - "null"
      - boolean
    doc: |
      reference sequences given on cmd line (as <reference_in>)
    inputBinding:
      position: 2
      prefix: '-c'

  large_index:
    type:
      - "null"
      - boolean
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides
    inputBinding:
      position: 3
      prefix: '--large-index'

  noauto:
    type:
      - "null"
      - boolean
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting
    inputBinding:
      position: 4
      prefix: '--noauto'

  packed:
    type:
      - "null"
      - boolean
    doc: |
      use packed strings internally; slower, less memory
    inputBinding:
      position: 5
      prefix: '--packed'

  bmax:
    type:
      - "null"
      - int
    doc: |
      max bucket sz for blockwise suffix-array builder
    inputBinding:
      position: 6
      prefix: '--bmax'

  bmaxdivn:
    type:
      - "null"
      - int
    doc: |
      max bucket sz as divisor of ref len (default: 4)
    inputBinding:
      position: 7
      prefix: '--bmaxdivn'

  dcv:
    type:
      - "null"
      - int
    doc: |
      diff-cover period for blockwise (default: 1024)
    inputBinding:
      position: 8
      prefix: '--dcv'

  nodc:
    type:
      - "null"
      - boolean
    doc: |
      disable diff-cover (algorithm becomes quadratic)
    inputBinding:
      position: 9
      prefix: '--nodc'

  noref:
    type:
      - "null"
      - boolean
    doc: |
      don't build .3/.4 index files
    inputBinding:
      position: 10
      prefix: '--noref'

  justref:
    type:
      - "null"
      - boolean
    doc: |
      just build .3/.4 index files
    inputBinding:
      position: 11
      prefix: '--justref'

  offrate:
    type:
      - "null"
      - int
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)
    inputBinding:
      position: 12
      prefix: '--offrate'

  ftabchars:
    type:
      - "null"
      - int
    doc: |
      # of chars consumed in initial lookup (default: 10)
    inputBinding:
      position: 13
      prefix: '--ftabchars'

  threads:
    type:
      - "null"
      - int
    doc: |
      # of threads
    inputBinding:
      position: 14
      prefix: '--threads'

  seed:
    type:
      - "null"
      - int
    doc: |
      seed for random number generator
    inputBinding:
      position: 15
      prefix: '--seed'

  quiet:
    type:
      - "null"
      - boolean
    doc: |
      verbose output (for debugging)
    inputBinding:
      position: 15
      prefix: '--quiet'

outputs:

  indices:
    type: File
    outputBinding:
      glob: $(inputs.bt2_index_base + ".1.bt2*")
    secondaryFiles: |
      ${
        var ext = self.location.split('/').slice(-1)[0].split('.').slice(-1)[0];
        var basename = self.location.split("/").slice(-1)[0].split(".").slice(0, -2).join (".");
        var dirname = self.location.split("/").slice(0,-1).join("/");
        return [{"class": "File", "path": dirname + "/" + basename + ".2." + ext},
                {"class": "File", "path": dirname + "/" + basename + ".3." + ext},
                {"class": "File", "path": dirname + "/" + basename + ".4." + ext},
                {"class": "File", "path": dirname + "/" + basename + ".rev.1." + ext},
                {"class": "File", "path": dirname + "/" + basename + ".rev.2." + ext}
        ]
      }

  output_log:
    type: File
    outputBinding:
      glob: $(inputs.bt2_index_base + ".log")

baseCommand:
  - bowtie2-build

arguments:
  - valueFrom: $('2> ' + inputs.bt2_index_base + '.log')
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bowtie2_build"
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool runs bowtie2-build to generate indices from input FASTA files