#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand: [samtools, stats]
requirements:
- class: DockerRequirement
  dockerPull: biocontainers/samtools:v1.7.0_cv3

inputs:
  input_file:
    type: File
    format: 
      - edam:format_2572  # BAM
      - edam:format_2573  # SAM
      - edam:format_3462  # CRAM
    inputBinding:
      position: 100
  coverage:
    type:
      - 'null'
      - type: record
        name: coverage_parameters
        fields:
          min_cov:
            type: int
          max_cov:
            type: int
          step_cov:
            type: int
    inputBinding:
      prefix: --coverage
    doc: "Set coverage distribution to the specified range (MIN, MAX, STEP all given as integers) [1,1000,1]"
  remove_dups:
    type: boolean?
    doc: "Exclude from statistics reads marked as duplicates"
    inputBinding:
      prefix: --remove_dups

  required_flag:
    type:
      - string
      - int
      - "null"
    default: 0
    doc: " STR|INT Required flag, 0 for unset. See also `samtools flags` [0] "
    inputBinding:
      prefix: -f

  filtering_flag:
    type:
      - string
      - int
      - "null"
    default: 0
    doc: "STR|INT Filtering flag, 0 for unset. See also `samtools flags` [0] "
    inputBinding:
      prefix: -F

  GC_depth:
    type: float?
    doc: "the size of GC-depth bins (decreasing bin size increases memory requirement) [2e4] "
    inputBinding:
      prefix: --GC-depth

  max_insert_size:
    type: int?
    doc: "Maximum insert size [8000]"
    inputBinding:
      prefix: -i

  listed_group:
    type: string?
    doc: "Include only listed read group or sample name [] "
    inputBinding:
      prefix: --id

  read_length:
    type: int?
    doc: "Include in the statistics only reads with the given read length [-1]"
    inputBinding:
      prefix: -l

  most_inserts:
    type: float?
    doc: "Report only the main part of inserts [0.99] "
    inputBinding:
      prefix: -m

  split_prefix:
    type: string?
    doc:  "A path or string prefix to prepend to filenames output when creating categorised statistics files with -S/--split. [input filename]"
    inputBinding:
      prefix: -P

  trim_quality:
    type: int?
    doc: "The BWA trimming parameter [0] "
    inputBinding:
      prefix: -q

  ref_seq:
    type: File?
    doc: "Reference sequence (required for GC-depth and mismatches-per-cycle calculation). [] "
    inputBinding:
      prefix: -r

  split:
    type: string?
    doc: "In addition to the complete statistics, also output categorised statistics based on the tagged field TAG (e.g., use --split RG to split into read groups).    Categorised statistics are written to files named <prefix>_<value>.bamstat, where prefix is as given by --split-prefix (or the input filename by default) and value has been encountered as the specified tagged field's value in one or more alignment records. "
    inputBinding:
      prefix: --split

  target_regions:
    type: File?
    doc: "Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive. []"
    inputBinding:
      prefix: --target-regions
  sparse:
    type: boolean?
    doc: "Suppress outputting IS rows where there are no insertions."
    inputBinding:
      prefix: --sparse
  remove_overlaps:
    type: boolean?
    doc: "Remove overlaps of paired-end reads from coverage and base count computations. "
    inputBinding:
      prefix: --remove-overlaps
  cov_threshold:
    type: int?
    doc: "Only bases with coverage above this value will be included in the target percentage computation [0] "
    inputBinding:
      prefix: -g

arguments:
  - prefix: --threads
    valueFrom: $(runtime.cores)  
    
# -X
#     If this option is set, it will allows user to specify customized index file location(s) if the data folder does not contain any index file. Example usage: samtools stats [options] -X /data_folder/data.bam /index_folder/data.bai chrM:1-10

outputs:
  stats:
    type: File
    outputBinding:
      glob: $(inputs.input_file.nameroot).stats.txt
stdout: $(inputs.input_file.nameroot).stats.txt

$namespaces: { edam: http://edamontology.org/, iana: https://www.iana.org/assignments/media-types/ }
$schemas:
  - 'http://edamontology.org/EDAM_1.18.owl'
