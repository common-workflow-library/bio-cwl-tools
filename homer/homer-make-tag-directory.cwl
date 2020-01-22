#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {"class": "Directory",
         "basename": "default",
         "listing": [inputs.bam_file],
         "writable": true}
      ]
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2

inputs:

  bam_file:
    type: File
    doc: "Alignment file, BAM"

  fragment_size:
    type:
      - "null"
      - int
      - string
    inputBinding:
      position: 5
      prefix: "-fragLength"
    doc: |
      Set fragment size.
      By default is estimated as if it was single end ChIP-Seq experiment.
      Possible values:
        "#" - int value to be used as fragment size
        "given" - use read lengths
        "pe" - calculate from paired end read coordinates

  total_reads:
    type:
      - "null"
      - int
      - string
    inputBinding:
      position: 6
      prefix: "-totalReads"
    doc: |
      Set total reads number for downstream normalization.
      Default: autocalculated, equal to uniquely mapped reads number
      Possible values:
        "#" - int value to be used as total reads number
        "all" - autocalculated, equal to uniquely + multi mapped reads number

  min_length:
    type: int?
    inputBinding:
      position: 7
      prefix: "-minlen"
    doc: |
      Discard reads smaller then

  max_length:
    type: int?
    inputBinding:
      position: 8
      prefix: "-maxlen"
    doc: |
      Discard reads bigger then

outputs:

  output_tag_folder:
    type: Directory
    outputBinding:
      glob: $(inputs.bam_file.basename.split('.')[0])
    doc: "Tag directory"

baseCommand: ["makeTagDirectory"]
arguments:
  - valueFrom: $(inputs.bam_file.basename.split('.')[0])
  - valueFrom: $("default/" + inputs.bam_file.basename)

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "homer-make-tag-directory"
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
  Tool runs makeTagDirectory that basically parses through the alignment file and splits the tags into separate
  files based on the chromosomes.

  Multiple alignment files are not supported. Alignment file's format is restricted to be only BAM.

  Output is placed in a folder with the name derived from the input BAM file's basename.

  Skipped arguments:

    Rely on the default value:
      -format           - format will be autodetected
      -precision        - the default value is used

    Not required general functionality:
      -d
      -single
      -force5th
      -t
      -flip
      -tbp

    Not required GC-bias options:
      -genome
      -checkGC
      -normGC
      -normFixedOligo
      -minNormRatio
      -maxNormRatio
      -iterNorm
      -filterReads

    Not required HiC options:
      -removePEbg
      -restrictionSite
      -removeSpikes
      -bowtiePE
      -directional