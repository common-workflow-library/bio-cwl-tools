#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7

doc: |
  sorting a bedgraph file by genomic coordinates
  
baseCommand: ["bash", "-c"]
arguments:
  - LC_COLLATE=C sort -k1,1 -k2,2n $(inputs.bedgraph.path)
stdout: |
  ${
    if( inputs.output_name == null ){
      return inputs.bedgraph.basename;
    }
    else{
      return inputs.output_name;
    }
  }

inputs:
  bedgraph:
    type: File
  output_name:
    doc: optional, if not specified, output file will be named as input file
    type: string?

outputs:
  bedgraph_sorted:
    type: stdout
    
