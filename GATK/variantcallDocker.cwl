#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: refIndex
  type: File
- id: refDict
  type: File
- id: referenceGenome
  type: File
- id: gatk_docker
  type: string
- id: gatk_path
  type: string
- id: name
  type: string
- id: AlignmentDocker_leftFastq
  type: File
- id: AlignmentDocker_rightFastq
  type: File

outputs:
- id: AlignmentDocker_rawSAM
  type: File
  outputSource: '#AlignmentDocker/rawSAM'
- id: AddOrReplaceReadGroupsDocker_rawBAM
  type: File
  outputSource: '#AddOrReplaceReadGroupsDocker/rawBAM'
- id: SortSamDocker_rawBAM
  type: File
  outputSource: '#SortSamDocker/rawBAM'
- id: ReferenceSeqIndexDocker_refIndex
  type: File
  outputSource: '#ReferenceSeqIndexDocker/refIndex'
- id: ReferenceSeqDictionaryDocker_refDict
  type: File
  outputSource: '#ReferenceSeqDictionaryDocker/refDict'
- id: MarkDuplicatesDocker_rawBAM
  type: File
  outputSource: '#MarkDuplicatesDocker/rawBAM'
- id: SplitNCigarReadsDocker_rawBAM
  type: File
  outputSource: '#SplitNCigarReadsDocker/rawBAM'
- id: HaplotypeCallerDocker_rawVCF
  type: File
  outputSource: '#HaplotypeCallerDocker/rawVCF'
- id: VariantFilterationDocker_rawVCF
  type: File
  outputSource: '#VariantFilterationDocker/rawVCF'
- id: SelectSNPsDocker_rawVCF
  type: File
  outputSource: '#SelectSNPsDocker/rawVCF'
- id: SelectINDELsDocker_rawVCF
  type: File
  outputSource: '#SelectINDELsDocker/rawVCF'

steps:
- id: AlignmentDocker
  in:
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  - id: index
    source: name
  - id: leftFastq
    source: AlignmentDocker_leftFastq
  - id: rightFastq
    source: AlignmentDocker_rightFastq
  run: AlignmentDocker.cwl
  out:
  - id: rawSAM
- id: AddOrReplaceReadGroupsDocker
  in:
  - id: inputSAM
    source: '#AlignmentDocker/rawSAM'
  - id: sampleName
    source: name
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  run: AddOrReplaceReadGroupsDocker.cwl
  out:
  - id: rawBAM
- id: SortSamDocker
  in:
  - id: inputBAM
    source: '#AddOrReplaceReadGroupsDocker/rawBAM'
  - id: sampleName
    source: name
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  run: SortSamDocker.cwl
  out:
  - id: rawBAM
- id: ReferenceSeqIndexDocker
  in:
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: ReferenceSeqIndexDocker.cwl
  out:
  - id: refIndex
- id: ReferenceSeqDictionaryDocker
  in:
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: ReferenceSeqDictionaryDocker.cwl
  out:
  - id: refDict
- id: MarkDuplicatesDocker
  in:
  - id: inputBAM
    source: '#SortSamDocker/rawBAM'
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: sampleName
    source: name
  run: MarkDuplicatesDocker.cwl
  out:
  - id: rawBAM
- id: SplitNCigarReadsDocker
  in:
  - id: inputBAM
    source: '#MarkDuplicatesDocker/rawBAM'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: SplitNCigarReadsDocker.cwl
  out:
  - id: rawBAM
- id: HaplotypeCallerDocker
  in:
  - id: inputBAM
    source: '#SplitNCigarReadsDocker/rawBAM'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: HaplotypeCallerDocker.cwl
  out:
  - id: rawVCF
- id: VariantFilterationDocker
  in:
  - id: mutantVCF
    source: '#HaplotypeCallerDocker/rawVCF'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: VariantFilterationDocker.cwl
  out:
  - id: rawVCF
- id: SelectSNPsDocker
  in:
  - id: mutantVCF
    source: '#VariantFilterationDocker/rawVCF'
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: sampleName
    source: name
  run: SelectSNPsDocker.cwl
  out:
  - id: rawVCF
- id: SelectINDELsDocker
  in:
  - id: mutantVCF
    source: '#VariantFilterationDocker/rawVCF'
  - id: docker
    source: gatk_docker
  - id: gatk_path
    source: gatk_path
  - id: ReferenceGenome
    source: referenceGenome
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: sampleName
    source: name
  run: SelectINDELsDocker.cwl
  out:
  - id: rawVCF
id: GATK4Docker
