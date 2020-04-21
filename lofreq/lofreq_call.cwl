#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: "LoFreq is a fast and sensitive variant-caller for inferring SNVs and indels
from next-generation sequencing data."
label: "LoFreq Call Variants"
hints:
  RessourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/lofreq:2.1.4--py27hc3dfafe_1

requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference_index)
      - $(inputs.reference_fasta)
      - $(inputs.reads_index)
      - $(inputs.reads_align)

baseCommand: [lofreq, call-parallel]

arguments:
  - prefix: --out
    valueFrom: $(inputs.reads_align.nameroot)_variant.vcf
    position: 99

inputs:
  - id: threads
    type: int?
    default: 1
    inputBinding:
      prefix: --pp-threads
      position: 1
  - id: reference_index
    type: File

  - id: reference_fasta
    doc: 'fasta'
    type: File
    format: edam:format_1930  # FASTA
    inputBinding:
      prefix: -f
      position: 1000
      valueFrom: $(self.basename)

  - id: call-indels
    type: boolean?
    inputBinding:
      prefix: --call-indels
      position: 3
    doc: "Enable indel calls (note: preprocess your file to include indel alignment qualities!)"

  - id: only-indels
    type: boolean?
    inputBinding:
      prefix: --only-indels
      position: 4
    doc: "Only call indels; no SNVs"

  - id: bed
    label: regions_from_bed
    type: File?
    doc: 'List of positions (chr pos) or regions (BED)'
    inputBinding:
      prefix: --bed

  - id: region
    type: string?
    doc: 'Limit calls to this region (chrom:start-end)'
    inputBinding:
      prefix: --region

  - id: min_bq
    label: min_base_quality
    type: int?
    default: 6
    inputBinding:
      prefix: --min-bq
    doc: 'Skip any base with baseQ smaller than INT [6]'

  - id: min-alt-bq
    label: min_alterne_base_quality
    type: int?
    default: 6
    inputBinding:
      prefix: --min-alt-bq
    doc: 'Skip alternate bases with baseQ smaller than INT [6]'

  - id: def-alt-bq
    label: def_alt_base_quality
    type: int?
    default: 0
    inputBinding:
      prefix: --def-alt-bq
    doc: 'Overwrite baseQs of alternate bases (that passed bq filter) with this value (-1: use median ref-bq; 0: keep) [0]'

  - id: min-jq
    label: min_joinedq
    type: int?
    default: 0
    inputBinding:
      prefix:  --min-jq
    doc: 'Skip any base with joinedQ smaller than INT [0]'

  - id: min-alt-jq
    label: min_alt_joinedq
    type: int?
    default: 0
    inputBinding:
      prefix: --min-alt-jq
    doc: "Skip alternate bases with joinedQ smaller than INT [0]"

  - id: def-alt-jq
    label: def_alt_joinedq
    type: int?
    default: 0
    inputBinding:
      prefix: --def-alt-jq
    doc: "Overwrite joinedQs of alternate bases (that passed jq filter) with this value (-1: use median ref-bq; 0: keep) [0]"

  - id: no-baq
    label: disable_base_alignment_quality
    type: boolean?
    inputBinding:
      prefix: --no-baq
    doc: 'Disable use of base-alignment quality (BAQ)'

  - id: no-idaq
    label: disable_indel_alignment_quality
    type: boolean?
    inputBinding:
      prefix: --no-idaq
    doc: "Don't use IDAQ values (NOT recommended under ANY circumstances other than debugging)"

  - id: del-baq
    label: delete_base_alignment_quality
    type: boolean?
    inputBinding:
      prefix: --del-baq
    doc: "Delete pre-existing BAQ values, i.e. compute even if already present in BAM"

  - id: no_ext_base_alignment_quality
    type: boolean?
    inputBinding:
      prefix: --no-ext-baq
    doc: "Use 'normal' BAQ (samtools default) instead of extended BAQ (both computed on the fly if not already present in lb tag)"

  - id: min_mq
    label: min_mapping_quality
    type: int?
    default: 0
    inputBinding:
      prefix: --min-mq
    doc: "Skip reads with mapping quality smaller than INT [0]"

  - id: max_mapping_quality
    type: int?
    default: 255
    inputBinding:
      prefix: --max-mq
    doc: "Cap mapping quality at INT [255]"

  - id: no_mapping_quality
    type: boolean?
    inputBinding:
      prefix: --no-mq
    doc: "Don't merge mapping quality in LoFreq's model"

  - id: enable_source_qual
    type: boolean?
    inputBinding:
      prefix: --src-qual
    doc: 'Enable computation of source quality'

  - id: ignore_vcf
    type: File[]?
    inputBinding:
      prefix: --ign-vcf
    doc: "Ignore variants in this vcf file for source quality computation. Multiple files can be given separated by commas"

  - id: replace_non_match
    type: int?
    default: -1
    inputBinding:
      prefix: --def-nm-q
    doc: 'If >= 0, then replace non-match base qualities with this default value [-1]'

  - id: pvalue_cutoff
    type: float?
    default: 0.01
    inputBinding:
      prefix: --sig
    doc: "P-Value cutoff / significance level [0.010000]"

  - id: bonferroni
    type: string?
    default: 'dynamic'
    inputBinding:
      prefix: --bonf
    doc: "Bonferroni factor. 'dynamic' (increase per actually performed test) or INT ['dynamic']"

  - id: min_cov
    type: int?
    default: 10
    inputBinding:
      prefix: --min-cov
      position: 2
    doc: "Test only positions having at least this coverage [1] (note: without --no-default-filter default filters (incl. coverage) kick in after predictions are done)"

  - id: max_depth_cov
    type: int?
    default: 1000000
    inputBinding:
      prefix: --max-depth
    doc: "Cap coverage at this depth [1000000]"

  - id: illumina_1_3
    type: boolean?
    inputBinding:
      prefix: --illumina-1.3
    doc: "Assume the quality is Illumina-1.3-1.7/ASCII+64 encoded"

  - id: use-orphan
    type: boolean?
    inputBinding:
      prefix: --use-orphan
    doc: "Count anomalous read pairs (i.e. where mate is not aligned properly)"

  - id: no-default-filter
    type: boolean?
    doc: "Don't run default lofreq filter automatically after calling variants"
    inputBinding:
      prefix: --no-default-filter

    # --plp-summary-only      No variant calling. Just output pileup summary per column
    # --force-overwrite       Overwrite any existing output
    # --verbose               Be verbose
    # --debug                 Enable debugging

  - id: reads_align
    doc: 'bam'
    type: File
    format: edam:format_2572  # BAM
    inputBinding:
      position: 1001
      valueFrom: $(self.basename)

  - id: reads_index
    doc: 'bai'
    type: File

outputs:
  - id: vcf
    type: File
    outputBinding:
      glob: "*.vcf"

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
