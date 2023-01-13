#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: biowardrobe2/rose:v0.0.2

inputs:
  binding_sites_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "GFF or BED file of binding sites used to make enhancers"

  bam_file:
    type: File
    inputBinding:
      position: 6
      prefix: "-r"
    secondaryFiles:
    - .bai
    doc: "Indexed bamfile to rank enhancer by"

  annotation_file:
    type: File
    inputBinding:
      position: 7
      prefix: "-g"
    doc: "TSV genome annotation file"

  stitch_distance:
    type: int
    inputBinding:
      position: 8
      prefix: "-s"
    doc: "Linking distance for stitching"

  tss_distance:
    type: int
    inputBinding:
      position: 9
      prefix: "-t"
    doc: "Distance from TSS to exclude. 0 = no TSS exclusion"

outputs:
  gff_directory:
    type: Directory
    outputBinding:
      glob: "gff"

  mapped_gff_directory:
    type: Directory
    outputBinding:
      glob: "mappedGFF"

  stitched_enhancer_region_map:
    type: File?
    outputBinding:
      glob: "*STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt"

  all_enhancers_table:
    type: File?
    outputBinding:
      glob: "*AllEnhancers.table.txt"

  super_enhancers_table:
    type: File?
    outputBinding:
      glob: "*SuperEnhancers.table.txt"

  enhancers_with_super_bed:
    type: File?
    outputBinding:
      glob: "*Enhancers_withSuper.bed"

  plot_points_pic:
    type: File?
    outputBinding:
      glob: "*Plot_points.png"

  gateway_enhancers_bed:
    type: File?
    outputBinding:
      glob: "*Gateway_Enhancers.bed"

  gateway_super_enhancers_bed:
    type: File?
    outputBinding:
      glob: "*Gateway_SuperEnhancers.bed"


baseCommand: ['ROSE_main', '-o', './']


$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

s:name: "rose"
s:license: http://www.apache.org/licenses/LICENSE-2.0

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
  Tool runs ROSE to get Super Enhancers regions
  -b and -c arguments are not supported
