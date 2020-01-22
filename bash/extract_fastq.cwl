#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Tool to decompress input FASTQ file
  Bash script's logic:
  - disable case sensitive glob check
  - check if root name of input file already include '.fastq' or '.fq' extension. If yes, set DEFAULT_EXT to ""
  - check file type, decompress if needed
  - return 1, if file type is not recognized
  This script also works of input file doesn't have any extension at all

requirements:
- class: ShellCommandRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3

inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      shopt -s nocaseglob

      FILE=$0
      DEFAULT_EXT=$1

      EXT_LIST=( ".fastq" ".fq" )

      BASENAME=$(basename "$FILE")
      ROOT_NAME="${BASENAME%.*}"

      for ITEM in $EXT_LIST; do
        if [[ $ROOT_NAME == *$ITEM ]]; then
          DEFAULT_EXT=""
        fi
      done

      T=`file -b "${FILE}" | awk '{print $1}'`
      case "${T}" in
        "bzip2"|"gzip"|"Zip")
          7z e -so "${FILE}" > "${ROOT_NAME}${DEFAULT_EXT}"
          ;;
        "ASCII")
          cp "${FILE}" "${ROOT_NAME}${DEFAULT_EXT}" || true
          ;;
        *)
          echo "Error: file type unknown"
          exit 1
      esac
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed FASTQ file

  compressed_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Compressed or uncompressed FASTQ file

  output_file_ext:
    type: string?
    inputBinding:
      position: 7
    default: ".fastq"
    doc: |
      Default extension for the extracted file

outputs:

  fastq_file:
    type: File
    outputBinding:
      glob: "*"

baseCommand: [bash, '-c']

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "extract_fastq"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bash/extract_fastq.cwl
s:codeRepository: https://github.com/common-workflow-library/bio-cwl-tools
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

s:about: |
  Tool to decompress input FASTQ file
