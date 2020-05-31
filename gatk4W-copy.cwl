#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow



requirements:

  StepInputExpressionRequirement: {}



doc: |

  Author: AMBARISH KUMAR er.ambarish@gmail.com & ambari73_sit@jnu.ac.in

  This is a proposed standard operating procedure for genomic variant detection using GATK4.

  It is hoped to be effective and useful for getting SARS-CoV-2 genome variants.

  

  It uses Illumina RNASEQ reads and genome sequence.



inputs:

  sars_cov_2_reference_genome:

    type: File

    


  sars_cov_2_reference_2bit_genome:

    type: File


  rnaseq_left_reads:

    type: File

    



  rnaseq_right_reads:

    type: File


  sampleName:

    type: string

   



steps:

  index_reference_genome_with_bowtie2:

    run: /home/ngsap2/tools/bowtie2/bowtie2_build.cwl

    in:

      reference_in: sars_cov_2_reference_genome

      bt2_index_base:

        valueFrom: "sars-cov-2"

    out: [ indices ]



  align_rnaseq_reads_to_genome:

    run: /home/ngsap2/tools/bowtie2/bowtie2_align.cwl

    in:

      indices_file: index_reference_genome_with_bowtie2/indices

      filelist: rnaseq_left_reads

      filelist_mates: rnaseq_right_reads

      output_filename:

        valueFrom: sars-cov-2.sam

    out: [ output ]



  index_reference_genome_with_samtools:

    run: /home/ngsap2/tools/samtools/samtools_faidx.cwl

    in:

      sequences: sars_cov_2_reference_genome

    out: [sequences_with_index]



  create_sequence_dictionary:

    run: /home/ngsap2/tools/picard/picard_CreateSequenceDictionary.cwl

    in:

      REFERENCE: index_reference_genome_with_samtools/sequences_with_index

    out: [ sequences_with_dictionary ]



  update_read_group:

    run: /home/ngsap2/tools/picard/picard_AddOrReplaceReadGroups.cwl

    in:

      INPUT: align_rnaseq_reads_to_genome/output

      OUTPUT:

        valueFrom: sars-cov-2-newreadgroups.bam

      RGID:

        valueFrom: "1"

      RGLB:

        valueFrom: 445_LIB

      RGPL:

        valueFrom: illumina

      RGSM:

        valueFrom: RNA

      RGPU:

        valueFrom: illumina

      SORT_ORDER:

        valueFrom: coordinate

    out: [ sequences_with_new_read_group ]

 

  sort_sam:
   
    run: /home/ngsap2/tools/gatk4/SortSamSparkNew.cwl

    in:
      
      inputBAM: update_read_group/sequences_with_new_read_group

      sampleName: sampleName

    out: [rawBAM] 

  
  mark_duplicates:

    run: /home/ngsap2/tools/gatk4/MarkDuplicatesSparkNew.cwl

    in:

      inputBAM: sort_sam/rawBAM
    
      sampleName: sampleName

    out: [ rawBAM ]



  split_alignments:

    run: /home/ngsap2/tools/GATK/GATK-SplitNCigarReads.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      reads: mark_duplicates/rawBAM

      output_filename:

        valueFrom: sars-cov-2-mutantsplit.bam

      # read_filter:  # Not available in GATK4

      #   valueFrom: ReassignOneMappingQuality  

    out: [ output ]



  index_split_alignments:

    run: /home/ngsap2/tools/samtools/samtools_index.cwl

    in:

      bam_sorted: split_alignments/output

    out: [ bam_sorted_indexed ]



  call_plausible_haplotypes_and_detect_variants:

    run: /home/ngsap2/tools/gatk4/HaplotypeCallerSparkNew.cwl

    in:

      RefDict: create_sequence_dictionary/sequences_with_dictionary
      
      Reference2bitGenome: sars_cov_2_reference_2bit_genome
      
      inputBAM: split_alignments/output
      
      RefIndex: index_reference_genome_with_samtools/sequences_with_index

      BAMindex: index_split_alignments/bam_sorted_indexed 

      sampleName: sampleName

      output_filename:

        valueFrom: sars-cov-2.mutant.vcf

    out: [rawVCF]



  filer_out_low_quality_variants:

    run: /home/ngsap2/tools/GATK/GATK-VariantFiltration.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      variant: call_plausible_haplotypes_and_detect_variants/rawVCF

      output_filename:

        valueFrom: sars-cov-2-mutantfilter.vcf

    out: [output]



  select_indels:

    run: /home/ngsap2/tools/GATK/GATK-SelectVariants.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      variant: filer_out_low_quality_variants/output

      select_type_to_include:

        valueFrom: INDEL

      output_filename:

        valueFrom: sars-cov-2-indel.vcf

    out: [ output ]



  select_snps:

    run: /home/ngsap2/tools/GATK/GATK-SelectVariants.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      variant: filer_out_low_quality_variants/output

      select_type_to_include:

        valueFrom: SNP

      output_filename:

        valueFrom: sars-cov-2-indel.vcf

    out: [ output ]

 

outputs:

  indels:

    type: File

    outputSource: select_indels/output

  snps:

    type: File

    outputSource: select_snps/output




