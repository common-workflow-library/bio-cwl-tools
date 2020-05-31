#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow



requirements:

  StepInputExpressionRequirement: {}



doc: |

  Author: AMBARISH KUMAR er.ambarish@gmail.com; ambari73_sit@jnu.ac.in

  This is a proposed standard operating procedure for genomic variant detection using VARSCAN.

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


  sample_name:

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



  sam_to_bam_conversion_using_samtools_view:

    run: /home/ngsap2/tools/samtools/samtools_view_sam2bam.cwl

    in: 

      sam: align_rnaseq_reads_to_genome/output

    out: [bam] 


  sort_alignment_files_using_samtools_sort:
   
    run: /home/ngsap2/tools/samtools/samtools_sort.cwl

    in: 

      bam_unsorted: sam_to_bam_conversion_using_samtools_view/bam

    out: [bam_sorted]


  
  index_bam_files_using_samtools_index:
    
    run: /home/ngsap2/tools/samtools/samtools_index.cwl

    in: 

      bam_sorted: sort_alignment_files_using_samtools_sort/bam_sorted

    out: [bam_sorted_indexed]


  mpileup_generation_using_samtools_mpileup:
   
    run: /home/ngsap2/tools/samtools/samtools_mpileup.cwl
    
    in:
      
      inputBAM: sort_alignment_files_using_samtools_sort/bam_sorted
      
      ReferenceGenome: sars_cov_2_reference_genome
      
      sampleName: sample_name 
      

    out: [rawMpileup]


  
  snp_calling_using_mpileup2snp:

    run: /home/ngsap2/tools/varscan/mpileup2snp.cwl

    in:

      inputMpileup: mpileup_generation_using_samtools_mpileup/rawMpileup

      sampleName: sample_name

    out: [snpVCF]



  indel_calling_using_mpileup2indel:

    run: /home/ngsap2/tools/varscan/mpileup2indel.cwl

    in:

      inputMpileup: mpileup_generation_using_samtools_mpileup/rawMpileup

      sampleName: sample_name

    out: [indelVCF]


outputs:

  snps:

    type: File

    outputSource: snp_calling_using_mpileup2snp/snpVCF

  indels:

    type: File

    outputSource: indel_calling_using_mpileup2indel/indelVCF 
   
$namespaces:

  edam: http://edamontology.org/

$schemas:

  - http://edamontology.org/EDAM_1.18.owl
