#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: ExpressionTool

# Inspired by https://github.com/galaxyproject/tools-iuc/blob/master/macros/read_group_macros.xml

requirements:
  InlineJavascriptRequirement:
    expressionLib:
     - |
         function clean(name) {
           return name.replace(/[^\w\-_\.]/, '_');
         };
         function read_group_name_default(input1, input2) {
           // input1: File
           // input2: File?
           var input_name1 = clean(input1.nameroot);
           if (input2 === undefined) {
             return input_name1;
           }
           var input_name2 = clean(input2.nameroot);
           var common_prefix = "";
           for (var index = 0; index < input_name1.length; index++) {
             if (input_name1.charAt(index) == input_name2.charAt(index)) {
              common_prefix += input_name1.charAt(index);
             }
           }
           if (common_prefix.length > 3) {
             return common_prefix;
           }
           return input_name1;
         };

inputs:
  input1: File
  input2: File?

# '@RG\tID:bwa-mem-fastq'

expression: |
  ${
    var rg_auto_name = read_group_name_default(inputs.input1, inputs.input2);
    var rg_id = rg_auto_name;
    var rg_string = "@RG\\tID:" + rg_id;
    return {"read_group_name": rg_string};
   }

outputs:
  read_group_name: string
