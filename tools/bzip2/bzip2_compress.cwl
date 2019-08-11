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
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bzip2_metadata.yaml

s:name: "bzip2_compress"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bzip2/bzip2_compress.cwl
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

doc: |
  Tool compresses `input_file` to `*.bz2`.
  Output file has the same basename, as input file, but with updated `.bz2` extension. `bzip2` exports compressed
  output file alognside the input file. To prevent tool from failing, `input_file` should be staged into output
  directory using `"writable": true`. Setting `writable: true` makes cwl-runner to make a copy of input file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode)

s:about: |
   usage: bzip2 [flags and input files in any order]

   -h --help           print this message
   -d --decompress     force decompression
   -z --compress       force compression
   -k --keep           keep (don't delete) input files
   -f --force          overwrite existing output files
   -t --test           test compressed file integrity
   -c --stdout         output to standard out
   -q --quiet          suppress noncritical error messages
   -v --verbose        be verbose (a 2nd -v gives more)
   -L --license        display software version & license
   -V --version        display software version & license
   -s --small          use less memory (at most 2500k)
   -1 .. -9            set block size to 100k .. 900k
   --fast              alias for -1
   --best              alias for -9

   If invoked as `bzip2', default action is to compress.
              as `bunzip2',  default action is to decompress.
              as `bzcat', default action is to decompress to stdout.

   If no file names are given, bzip2 compresses or decompresses
   from standard input to standard output.  You can combine
   short flags, so `-v -4' means the same as -v4 or -4v, &c.
