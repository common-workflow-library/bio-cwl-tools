#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Kallisto index
doc: |

  Docs: https://pachterlab.github.io/kallisto/


  Builds a kallisto index

  Usage: kallisto index [arguments] FASTA-files

  Required argument:
  -i, --index=STRING          Filename for the kallisto index to be constructed 

  Optional argument:
  -k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 63)
  -t, --threads=INT           Number of threads to use (default: 1)
  -d, --d-list=STRING         Path to a FASTA-file containing sequences to mask from quantification
      --make-unique           Replace repeated target names with unique names
      --aa                    Generate index from a FASTA-file containing amino acid sequences
      --distinguish           Generate index where sequences are distinguished by the sequence name
  -T, --tmp=STRING            Temporary directory (default: tmp)
  -m, --min-size=INT          Length of minimizers (default: automatically chosen)
  -e, --ec-max-size=INT       Maximum number of targets in an equivalence class (default: no maximum)

  This CWL was adapted from: https://github.com/common-workflow-library/bio-cwl-tools/commit/91c42fb809ce18eafe16155cca0abf362270c0fe

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1
  SoftwareRequirement:
    packages:
      - package: kallisto
        version: [ "0.51.1" ]
        specs:
          - https://identifiers.org/rrid/RRID:SCR_016582
          - https://identifiers.org/biotools/kallisto 

baseCommand: [kallisto, index]

inputs:
  inputFasta:
    type: File[]
    format: edam:format_1929 # FASTA
    inputBinding:
      position: 200

  IndexName:
    type: string
    inputBinding:
      prefix: "--index="
      separate: false
      valueFrom: $(self).kallistoIndex

#Optional arguments

  kmerSize:
    type: int?
    inputBinding:
      prefix: "--kmer-size="
      separate: false

  makeUnique:
    type: boolean?
    inputBinding:
      prefix: "--make-unique"

outputs:
  index:
    type: File
    outputBinding:
      glob: $(inputs.IndexName).kallistoIndex

$namespaces:
  edam: https://edamontology.org/
  s: https://schema.org/
$schemas:
  - https://edamontology.org/EDAM_1.25.owl
  - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: https://spdx.org/licenses/BSD-2-Clause
s:citation: https://dx.doi.org/10.1038/nbt.3519
s:codeRepository: https://github.com/pachterlab/kallisto
