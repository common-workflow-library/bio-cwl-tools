#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/ucsc-bedgraphtobigwig:377--ha8a8165_3
  SoftwareRequirement:
    packages:
      bedgraphtobigwig:
        specs:
          - https://bio.tools/bedgraphtobigwig
          - https://anaconda.org/bioconda/ucsc-bedgraphtobigwig

inputs:
  bedgraph_file:
    type: File
    format: edam:format_3583  # bedGraph
    inputBinding:
      position: 10
    doc: |
      Four column bedGraph file: <chrom> <start> <end> <value>

  chrom_length_file:
    type: File
    inputBinding:
      position: 11
    doc: |
      Two-column chromosome length file: <chromosome name> <size in bases>

  unc:
    type: boolean?
    inputBinding:
      position: 5
      prefix: "-unc"
    doc: |
      Disable compression

  items_per_slot:
    type: int?
    inputBinding:
      separate: false
      position: 6
      prefix: "-itemsPerSlot="
    doc: |
      Number of data points bundled at lowest level. Default 1024

  block_size:
    type: int?
    inputBinding:
      separate: false
      position: 7
      prefix: "-blockSize="
    doc: |
      Number of items to bundle in r-tree.  Default 256

  output_filename:
    type: string?
    inputBinding:
      position: 12
      valueFrom: |
        $( self == "" ? inputs.bedgraph_file.nameroot + ".bigWig" : self )
    default: ""
    doc: |
      If set, writes the output bigWig file to output_filename,
      otherwise generates filename from the nameroot of the bedgraph_file

outputs:
  bigwig_file:
    type: File
    format: edam:format_3006  # BigWig
    outputBinding:
      glob: |
        $( self == "" ? inputs.bedgraph_file.nameroot + ".bigWig" : self )

baseCommand: bedGraphToBigWig

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "ucsc-bedgraphtobigwig"
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
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Tool converts bedGraph to bigWig file.

  `default_output_filename` function returns filename for generated bigWig if `output_filename` is not provided.
  Default filename is generated on the base of `bedgraph_file` basename with the updated to `*.bigWig` extension.
