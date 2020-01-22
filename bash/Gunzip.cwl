#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: ubuntu:xenial

inputs:
  InputFile:
    type: File
    inputBinding:
      position: 1

baseCommand: [ "gunzip" ]

arguments: [ "-c" ]

outputs:
  unzipped_file:
    type: stdout

stdout: unzippedfile.stdout
# Source used from https://github.com/griffithlab/pmbio.org/blob/master/assets/CWL/gunzip.cwl