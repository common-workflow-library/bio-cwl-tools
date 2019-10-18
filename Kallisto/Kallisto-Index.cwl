cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0"
  InlineJavascriptRequirement: {}


    
  IndexName:
    type: string
    inputBinding:
      prefix: "--index="
      separate: false
      valueFrom: $(self + ".kl")

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

baseCommand: [kallisto, index]

outputs:

  index:
    type: File
    outputBinding:
      glob: $(inputs.IndexName)
  