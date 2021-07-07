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

    format: edam:format_1929  # FASTA



  rnaseq_left_reads:

    type: File

    format: edam:format_1930  # FASTQ



  rnaseq_right_reads:

    type: File

    format: edam:format_1930  # FASTQ
    
  indices_folder:
  
    type: Directory
    



steps:

  index_reference_genome_with_bowtie2:

    run: bowtie2_build.cwl

    in:

      reference_in: sars_cov_2_reference_genome

      bt2_index_base:

        valueFrom: "sars-cov-2"

    out: [ indices ]



  align_rnaseq_reads_to_genome:

    run: bowtie2_align.cwl

    in:

      indices_file: index_reference_genome_with_bowtie2/indices

      filelist: rnaseq_left_reads

      filelist_mates: rnaseq_right_reads
      
      indices_folder: indices_folder

      output_filename:

        valueFrom: sars-cov-2.sam

    out: [ output ]



  index_reference_genome_with_samtools:

    run: samtools_faidx.cwl

    in:

      sequences: sars_cov_2_reference_genome

    out: [sequences_with_index]



  create_sequence_dictionary:

    run: picard_CreateSequenceDictionary.cwl

    in:

      REFERENCE: sars_cov_2_reference_genome

    out: [ sequences_with_dictionary ]



  update_read_group:

    run: picard_AddOrReplaceReadGroups.cwl

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

 

  mark_duplicates:

    run: picard_markdup.cwl

    in:

      bam_sorted: update_read_group/sequences_with_new_read_group

    out: [ bam_duprem ]



  split_alignments:

    run: GATK-SplitNCigarReads.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      reads: mark_duplicates/bam_duprem

      output_filename:

        valueFrom: sars-cov-2-mutantsplit.bam

      # read_filter:  # Not available in GATK4

      #   valueFrom: ReassignOneMappingQuality  

    out: [ output ]



  index_split_alignments:

    run: samtools_index.cwl

    in:

      bam_sorted: split_alignments/output

    out: [ bam_sorted_indexed ]



  call_plausible_haplotypes_and_detect_variants:

    run: GATK-HaplotypeCaller.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      input: index_split_alignments/bam_sorted_indexed

      output_filename:

        valueFrom: sars-cov-2-mutant.vcf

    out: [ output ]



  filer_out_low_quality_variants:

    run: GATK-VariantFiltration.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      variant: call_plausible_haplotypes_and_detect_variants/output

      output_filename:

        valueFrom: sars-cov-2-mutantfilter.vcf

    out: [output]



  select_indels:

    run: GATK-SelectVariants.cwl

    in:

      reference: create_sequence_dictionary/sequences_with_dictionary

      variant: filer_out_low_quality_variants/output

      select_type_to_include:

        valueFrom: INDEL

      output_filename:

        valueFrom: sars-cov-2-indel.vcf

    out: [ output ]



  select_snps:

    run: GATK-SelectVariants.cwl

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



$namespaces:

  edam: http://edamontology.org/

$schemas:

  - http://edamontology.org/EDAM_1.18.owl
