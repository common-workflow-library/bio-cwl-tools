#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: ExpressionTool

doc: |
  Also consider https://www.commonwl.org/user_guide/misc/#rename-an-input-file
  or https://www.commonwl.org/user_guide/misc/#rename-an-output-file

requirements:
  InlineJavascriptRequirement: {}

inputs:
  srcfile:
    type: File
  newname:
    type: string

expression: |
  ${
     inputs.srcfile.basename = inputs.newname;
     return {"outfile": inputs.srcfile};
   }

outputs:
  outfile:
    type: File
