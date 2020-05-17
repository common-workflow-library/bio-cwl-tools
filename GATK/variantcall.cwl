#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: gatk
  type: File
- id: refIndex
  type: File
- id: refDict
  type: File
- id: referenceGenome
  type: File
- id: name
  type: string
- id: Alignment_leftFastq
  type: File
- id: Alignment_rightFastq
  type: File

outputs:
- id: Alignment_rawSAM
  type: File
  outputSource: '#Alignment/rawSAM'
- id: AddOrReplaceReadGroups_rawBAM
  type: File
  outputSource: '#AddOrReplaceReadGroups/rawBAM'
- id: SortSam_rawBAM
  type: File
  outputSource: '#SortSam/rawBAM'
- id: ReferenceSeqIndex_refIndex
  type: File
  outputSource: '#ReferenceSeqIndex/refIndex'
- id: ReferenceSeqDictionary_refDict
  type: File
  outputSource: '#ReferenceSeqDictionary/refDict'
- id: MarkDuplicates_rawBAM
  type: File
  outputSource: '#MarkDuplicates/rawBAM'
- id: SplitNCigarReads_rawBAM
  type: File
  outputSource: '#SplitNCigarReads/rawBAM'
- id: HaplotypeCaller_rawVCF
  type: File
  outputSource: '#HaplotypeCaller/rawVCF'
- id: VariantFilteration_rawVCF
  type: File
  outputSource: '#VariantFilteration/rawVCF'
- id: SelectSNPs_rawVCF
  type: File
  outputSource: '#SelectSNPs/rawVCF'
- id: SelectINDELs_rawVCF
  type: File
  outputSource: '#SelectINDELs/rawVCF'

steps:
- id: Alignment
  in:
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  - id: index
    source: name
  - id: leftFastq
    source: Alignment_leftFastq
  - id: rightFastq
    source: Alignment_rightFastq
  run: Alignment.cwl
  out:
  - id: rawSAM
- id: AddOrReplaceReadGroups
  in:
  - id: inputSAM
    source: '#Alignment/rawSAM'
  - id: GATK
    source: gatk
  - id: sampleName
    source: name
  run: AddOrReplaceReadGroups.cwl
  out:
  - id: rawBAM
- id: SortSam
  in:
  - id: inputBAM
    source: '#AddOrReplaceReadGroups/rawBAM'
  - id: GATK
    source: gatk
  - id: sampleName
    source: name
  run: SortSam.cwl
  out:
  - id: rawBAM
- id: ReferenceSeqIndex
  in:
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: ReferenceSeqIndex.cwl
  out:
  - id: refIndex
- id: ReferenceSeqDictionary
  in:
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: ReferenceSeqDictionary.cwl
  out:
  - id: refDict
- id: MarkDuplicates
  in:
  - id: inputBAM
    source: '#SortSam/rawBAM'
  - id: GATK
    source: gatk
  - id: sampleName
    source: name
  run: MarkDuplicates.cwl
  out:
  - id: rawBAM
- id: SplitNCigarReads
  in:
  - id: inputBAM
    source: '#MarkDuplicates/rawBAM'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: SplitNCigarReads.cwl
  out:
  - id: rawBAM
- id: HaplotypeCaller
  in:
  - id: inputBAM
    source: '#SplitNCigarReads/rawBAM'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: HaplotypeCaller.cwl
  out:
  - id: rawVCF
- id: VariantFilteration
  in:
  - id: mutantVCF
    source: '#HaplotypeCaller/rawVCF'
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: sampleName
    source: name
  run: VariantFilteration.cwl
  out:
  - id: rawVCF
- id: SelectSNPs
  in:
  - id: mutantVCF
    source: '#VariantFilteration/rawVCF'
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: sampleName
    source: name
  run: SelectSNPs.cwl
  out:
  - id: rawVCF
- id: SelectINDELs
  in:
  - id: mutantVCF
    source: '#VariantFilteration/rawVCF'
  - id: GATK
    source: gatk
  - id: ReferenceGenome
    source: referenceGenome
  - id: RefIndex
    source: refIndex
  - id: RefDict
    source: refDict
  - id: sampleName
    source: name
  run: SelectINDELs.cwl
  out:
  - id: rawVCF
id: GATK4
