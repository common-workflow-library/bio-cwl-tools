class: CommandLineTool
cwlVersion: v1.0
id: picard_SortSam

hints:
  - class: DockerRequirement
    dockerPull: broadinstitute/picard:latest
requirements:
  InlineJavascriptRequirement: {}

baseCommand: [java, -jar, /usr/picard/picard.jar]
arguments:
  - SortSam
#  - VERBOSY=ERROR

inputs:
  - id: inputFile
    type: File
    inputBinding:
      position: 1
      prefix: INPUT=
      separate: false

  - id: outFile_name
    type: string?
    inputBinding:
      position: 2
      prefix: OUTPUT=
      valueFrom: |
        ${ if(inputs.sort_order == "coordinate") { return (inputs.inputFile.nameroot)+".bam";} else { return (inputs.inputFile.nameroot)+".sam"; } }
      separate: false

  - id: sort_order
    type: string?
    default: coordinate
    doc: 'coordinate (bam) or queryname (sam)'
    inputBinding:
      position: 6
      prefix: SORT_ORDER=
      separate: false

  - id: validation_stringency
    type: string?
    default: LENIENT
    doc: 'LENIENT, SILENT or STRICT'
    inputBinding:
      position: 11
      prefix: VALIDATION_STRINGENCY=
      separate: false


outputs:
  - id: outFile
    type: File
    outputBinding:
      glob: '*.*am'
