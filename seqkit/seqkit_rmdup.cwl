#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
hints:
  SoftwareRequirement:
    packages:
      seqkit:
        version: [ 0.7.1 ]
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/seqkit:0.7.1--0"
  ResourceRequirement:
    coresMin: 8
    coresMax: 32
    ramMin: $(7 * 1024)
    outdirMin: |
      ${
        var sum = 0;
        for (var i = 0; i < inputs.readsFA.length; i++) {
          sum += inputs.readsFA[i].size;
        }
        return (sum/(1024*1024*1024)+1) + 20;
      }

inputs:
  reads: File
  ignore_case:
    type: boolean?
    inputBinding:
      prefix: --ignore-case

outputs:
  reads_dedup:
    type: File
    outputBinding:
      glob: $(inputs.reads.nameroot)_dedup.fasta
  dups:
    type: File?
    outputBinding:
      glob: dups.txt

baseCommand: [ seqkit, rmdup ]

arguments:
 - --by-seq
 - --threads=$(runtime.cores)
 - --dup-num-file,
 - dups.txt,
 - --out-file,
 - $(inputs.reads.nameroot)_dedup.fasta
 - $(inputs.reads.path)
