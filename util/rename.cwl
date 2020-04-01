#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: "true"

requirements:
  InlineJavascriptRequirement: {}

inputs:
  srcfile:
    type: File
  newname:
    type: string

outputs:
  outfile:
    type: File
    outputBinding:
      outputEval: |
        ${
          console.log(self);
          var file = inputs.srcfile;
          file.basename = inputs.newname;
          return file;
        }

