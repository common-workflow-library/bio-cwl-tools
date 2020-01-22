#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2

inputs:

  peak_file:
    type: File
    doc: "Homer generated peak file or BED"

  tag_folders:
    type:
      - Directory
      - Directory[]
    inputBinding:
      position: 7
      prefix: "-d"
    doc: "Tag folders from homer-make-tag-directory tool"

  flank_width:
    type: int?
    inputBinding:
      position: 8
      prefix: "-size"
    doc: |
      Size of 5' and 3' flanks, default: 5000

  flank_bin_size:
    type: int?
    inputBinding:
      position: 9
      prefix: "-bin"
    doc: |
      Bin size for 5' and 3' flanks, default: 100

  genebody_bin_count:
    type: int?
    inputBinding:
      position: 10
      prefix: "-gbin"
    doc: |
      Number of bins in gene body, default: 50

  min_peak_size:
    type: int?
    inputBinding:
      position: 11
      prefix: "-min"
    doc: |
      Minimum size of peak region to use, default: 3000

  max_peak_size:
    type: int?
    inputBinding:
      position: 12
      prefix: "-max"
    doc: |
      Maximum size of peak region to use, default: no max

  genebody_to_flank_ratio:
    type: int?
    inputBinding:
      position: 13
      prefix: "-gRatio"
    doc: |
      Ratio of gene region to flanks, default: 2

  norm_fpkm:
    type: boolean?
    inputBinding:
      position: 14
      prefix: "-fpkm"
    doc: |
      Normalize read counts to million reads or fragments per kilobase mapped

  norm_raw:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "-raw"
    doc: |
      Do not adjust the tag counts based on total tags sequenced.
      By default all tag counts will be normalized to norm_tag_count

  norm_tag_count:
    type: int?
    inputBinding:
      position: 16
      prefix: "-norm"
    doc: |
      Normalize tags to this tag count, default=1e7, 0=average tag count in all directories

  norm_fragment_size:
    type: int?
    inputBinding:
      position: 17
      prefix: "-normLength"
    doc: |
      Fragment length to normlize to for experiments with different lens. Default: 100bp

  strand:
    type: string?
    inputBinding:
      position: 18
      prefix: "-strand"
    doc: |
      Count tags on specific strands relative to peak. Default: both
      Possible values: +|-

  threads:
    type: int?
    inputBinding:
      position: 19
      prefix: "-cpu"
    doc: |
      Set the number of threads. Default: 1

  histogram_filename:
    type: string
    doc: "Output histogram's filename"

outputs:

  histogram_file:
    type: stdout
    doc: "Output histogram file"

stdout: ${return inputs.histogram_filename;}

baseCommand: ["makeMetaGeneProfile.pl"]
arguments:
  - valueFrom: $(inputs.peak_file)
    position: 5
  - valueFrom: $("none")
    position: 6

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "homer-make-metagene-profile"
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
  Tool runs makeMetaGeneProfile.pl which automates the creation of 'meta-gene' profiles.