#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

doc: |
  Generates coverage tracks in bedgraph or bigiwig format from BAM file.
  Normalization by spike-in reads is supported.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1
  SoftwareRequirement:
    packages:
      deeptools:
        specs: [ "http://identifiers.org/biotools/deeptools" ]
        version: [ "3.1.1" ]

baseCommand: ["bamCoverage"]
arguments:  
  - valueFrom: $(inputs.bam.nameroot).bigwig
    prefix: --outFileName
    position: 10
  - valueFrom: |
        ${ 
          if( inputs.spike_in_count == null ){
            return inputs.normalizeUsing
          }
          else{
            return null 
          }
        }
    prefix: --normalizeUsing
    position: 10
  
inputs:
  bam:
    doc: bam file as input; needs bai index file in the same directory
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bam
  effective_genome_size:
    doc: |
      the effectively mappable genome size, 
      see: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    type: long
    inputBinding:
        position: 10
        prefix: --effectiveGenomeSize
  bin_size:
    type: int
    default: 10
    inputBinding:
      prefix: --binSize
      position: 10
  ignoreForNormalization:
    type:
      type: array
      items: string
    default: ["chrX", "chrY", "chrM"]
    inputBinding:
      prefix: --ignoreForNormalization
      position: 10
  normalizeUsing:
    type: string
    default: "RPGC"
  extendReads:
    doc: For single end reads, the fragment size should be specified here.
    type: int?
    inputBinding:
      prefix: --extendReads
      position: 10
  spike_in_count:
    doc: |
      Number of reads aligned to the spike in reference, optional.
      If specified, coverage will be multiplied by 1/spike_in_count and
      normalizeUsing will be ignored.
    type: long?
    inputBinding:
      position: 10
      prefix: --scaleFactor
      valueFrom: |
        ${ 
          if( self == null ){
            return null
          }
          else{
            return (1.0 / parseFloat(self)) 
          }
        }
  outFileFormat:
    doc: bigiwg or bedgraph
    type: string
    default: bigwig
    inputBinding:
      position: 10
      prefix: --outFileFormat

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot).bigwig
    
