#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Calls nucleosome positions based on ATAC-seq reads.

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/nucleoatac:0.3.4

baseCommand: ["nucleoatac", "run"]
stderr: $(inputs.output_basename).stderr
stdout: $(inputs.output_basename).stout

inputs:
  bam:
    doc: aligned and filtered reads, not shifted
    type: File
    inputBinding:
      prefix: --bam
      position: 1
  bed:
    doc: regions of open chromatin (ATAC peaks)
    type: File
    inputBinding:
      prefix: --bed
      position: 1
  fasta:
    doc: reference genome in fasta format + samtools faidx index
    type: File
    secondaryFiles: [ ".fai" ]
    inputBinding:
      prefix: --fasta
      position: 1
  output_basename:
    type: string
    inputBinding:
      prefix: --out
      position: 1

outputs:
  nucl_occ_tracks:
    type: File?
    outputBinding:
      glob:  "*.occ.bedgraph.gz"
  nucl_occ_lower_bound_tracks:
    type: File?
    outputBinding:
      glob:  "*.occ.lower_bound.bedgraph.gz"
  nucl_occ_upper_bound_tracks:
    type: File?
    outputBinding:
      glob:  "*.occ.upper_bound.bedgraph.gz"
  nucl_dist_txt:
    type: File?
    outputBinding:
      glob:  "*.nuc_dist.txt"
  nucl_dist_plot:
    type: File?
    outputBinding:
      glob:  "*.nuc_dist.eps"
  fragsize_in_peaks_txt:
    type: File?
    outputBinding:
      glob:  "*.fragmentsizes.txt"
  nucl_occ_fit_txt:
    type: File?
    outputBinding:
      glob:  "*.occ_fit.txt"
  nucl_occ_fit_plot:
    type: File?
    outputBinding:
      glob:  "*.occ_fit.eps"
  nucl_occ_peaks_bed:
    type: File?
    outputBinding:
      glob:  "*.occpeaks.bed.gz"
  nucl_vplot_data:
    type: File?
    outputBinding:
      glob:  "*.VMat"
  nucl_pos_bed:
    type: File?
    outputBinding:
      glob:  "*.nucpos.bed.gz"
  nucl_pos_redundant_bed:
    type: File?
    outputBinding:
      glob:  "*.nucpos.redundant.bed.gz"
  nucl_norm_crosscor_tracks:
    type: File?
    outputBinding:
      glob:  "*.nucleoatac_signal.bedgraph.gz"
  nucl_norm_smooth_crosscor_tracks:
    type: File?
    outputBinding:
      glob:  "*.nucleoatac_signal.smooth.bedgraph.gz"
  combined_nucl_pos_bed:
    type: File?
    outputBinding:
      glob:  "*.nucmap_combined.bed.gz"
  nfr_pos_bed:
    type: File?
    outputBinding:
      glob:  "*.nfrpos.bed.gz"
  nucleoatac_stderr:
    type: stderr
  nucleoatac_stdout:
    type: stdout
