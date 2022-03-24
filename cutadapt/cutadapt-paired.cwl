#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
   dockerPull: quay.io/biocontainers/cutadapt:3.7--py39hbf8eff0_1

baseCommand: cutadapt

inputs:
  reads_1:
    type: File
    inputBinding: {}
  reads_2:
    type: File
    inputBinding: {}
  minimum_length:
    type: int?
    inputBinding:
      prefix: --minimum-length

  quality_cutoff:
    type: int?
    inputBinding:
      prefix: --quality-cutoff

stdout: report.txt

arguments:
  - prefix: --output
    valueFrom: $(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)
  - prefix: --paired-output
    valueFrom: $(inputs.reads_2.basename).trimmed$(inputs.reads_2.nameext)

outputs:
  report: stdout
  trimmed_reads_1:
    type: File
    format: $(inputs.reads_1.format)
    outputBinding:
      glob: $(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)
  trimmed_reads_2:
    type: File
    format: $(inputs.reads_2.format)
    outputBinding:
      glob: $(inputs.reads_2.basename).trimmed$(inputs.reads_2.nameext)


