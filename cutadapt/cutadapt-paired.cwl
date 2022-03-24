#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
   dockerPull: quay.io/biocontainers/cutadapt:3.7--py39hbf8eff0_1

baseCommand: cutadapt

inputs:
  reads_1: File
  reads_2: File
  minimum_length:
    type: int?
    inputBinding:
      prefix: --minimum_length

  quality_cutoff:
    type: int?
    inputBinding:
      prefix: --quality-cutoff

stdout: report.txt

arguments:
  - prefix: --output
    valueFrom: $(inputs.reads_1.basename).trimmed.$(inputs.reads_1.nameext)
  - prefix: --paired-output
    valueFrom: $(inputs.reads_2.basename).trimmed.$(inputs.reads_2.nameext)

outputs:
  report: stdout
  trimmed_reads_1:
    type: File
    outputBinding:
      glob: $(inputs.reads_1.basename).trimmed.$(inputs.reads_1.nameext)
  trimmed_reads_2:
    type: File
    outputBinding:
      glob: $(inputs.reads_2.basename).trimmed.$(inputs.reads_2.nameext)


