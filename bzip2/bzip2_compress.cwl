#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.input_file,
                  "entryname": inputs.input_file.basename,
                  "writable": true
                }
              ]
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2

inputs:

  input_file:
    type:
      - File
    inputBinding:
      position: 1
    doc: |
      File to be compressed

  fast:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 2
    doc: |
      Set block size to 100k

  best:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 3
    doc: |
      Set block size to 900k

outputs:

  output_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return inputs.input_file.basename + '.bz2';
        }

baseCommand: [bzip2]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bzip2_compress"
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
  Tool compresses `input_file` to `*.bz2`.
  Output file has the same basename, as input file, but with updated `.bz2` extension. `bzip2` exports compressed
  output file alognside the input file. To prevent tool from failing, `input_file` should be staged into output
  directory using `"writable": true`. Setting `writable: true` makes cwl-runner to make a copy of input file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode)