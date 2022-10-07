#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/rgt:v0.0.1


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      add_rep () {
          echo "$1" >> "$3"
          while IFS=',' read -ra REP; do
              for i in "${REP[@]}"; do
                  echo $i >> "$3"
              done
          done <<< "$2"
      }
      > config.txt
      add_rep "#rep1" "$0" config.txt
      add_rep "#rep2" "$1" config.txt
      echo "#chrom_sizes" >> config.txt
      echo "$2" >> config.txt
      rgt-THOR config.txt "${@:3}"
    inputBinding:
      position: 5
    doc: "Bash function to generate config.txt and run rgt-THOR"

  bambai_pair_cond_1:
    type: File[]
    inputBinding:
      position: 6
      itemSeparator: ","
    secondaryFiles:
    - .bai
    doc: "Alignment and index files for the first biological condition"

  bambai_pair_cond_2:
    type: File[]
    inputBinding:
      position: 7
      itemSeparator: ","
    secondaryFiles:
    - .bai
    doc: "Alignment and index files for the second biological condition"

  chrom_length_file:
    type: File
    inputBinding:
      position: 8
    doc: "Chromosome length file"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--name"
      position: 9
    doc: "Experiment's name and prefix for all files that are created"

  merge_peaks:
    type: boolean?
    inputBinding:
      prefix: "--merge"
      position: 10
    doc: "Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data)"

  no_merge_bin:
    type: boolean?
    inputBinding:
      prefix: "--no-merge-bin"
      position: 11
    doc: "Merge the overlapping bin before filtering by p-value"

  housekeeping_genes_bed_file:
    type: File?
    inputBinding:
      prefix: "--housekeeping-genes"
      position: 12
    doc: "Define housekeeping genes (BED format) used for normalizing"

  deadzones_bed_file:
    type: File?
    inputBinding:
      prefix: "--deadzones"
      position: 14
    doc: "Define blacklisted genomic regions avoided for analysis"

  no_correction:
    type: boolean?
    inputBinding:
      prefix: "--no-correction"
      position: 15
    doc: "Do not use multipe test correction for p-values (Benjamini/Hochberg)"

  pvalue_cutoff:
    type: float?
    inputBinding:
      prefix: "--pvalue"
      position: 16
    doc: "P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: 0.1]"

  extension_size:
    type:
      - "null"
      - string
      - int[]
    inputBinding:
      prefix: "--exts"
      position: 17
      itemSeparator: ","
    doc: |
      Read's extension size for BAM files (comma separated list for each BAM file in config file).
      If option is not chosen, estimate extension sizes

  normalization_factor:
    type:
      - "null"
      - float[]
    inputBinding:
      prefix: "--factors-inputs"
      position: 18
      itemSeparator: ","
    doc: |
      Normalization factors for input-DNA (comma separated list for each BAM file in config file).
      If option is not chosen, estimate factors

  scaling_factor:
    type:
      - "null"
      - float[]
    inputBinding:
      prefix: "--scaling-factors"
      position: 19
      itemSeparator: ","
    doc: |
      Scaling factor for each BAM file (not control input-DNA) as comma separated list for each BAM file in config file.
      If option is not chosen, follow normalization strategy (TMM or HK approach) [default: none]

  regions_bed_file:
    type: File?
    inputBinding:
      prefix: "--regions"
      position: 20
    doc: "Define regions (BED format) to restrict the analysis, that is, where to train the HMM and search for DPs"

  bin_size:
    type: int?
    inputBinding:
      prefix: "--binsize"
      position: 21
    doc: "Size of underlying bins for creating the signal. [default: 100]"

  step_size:
    type: int?
    inputBinding:
      prefix: "--step"
      position: 22
    doc: "Stepsize with which the window consecutively slides across the genome to create the signal. [default: 50]"

  no_gc_content:
    type: boolean?
    inputBinding:
      prefix: "--no-gc-content"
      position: 23
    doc: "Do not normalize towards GC content"

  norm_regions_bed_file:
    type: File?
    inputBinding:
      prefix: "--norm-regions"
      position: 24
    doc: "Restrict normalization to particular regions"

  foldchange:
    type: float?
    inputBinding:
      prefix: "--foldchange"
      position: 25
    doc: "Fold change parameter to define training set (t_1, see paper). [default: 1.6]"

  threshold:
    type: int?
    inputBinding:
      prefix: "--threshold"
      position: 26
    doc: "Minimum signal support for differential peaks to define training set as percentage (t_2, see paper). [default: 95]"

  bin_count:
    type: int?
    inputBinding:
      prefix: "--size"
      position: 27
    doc: "Number of bins the HMM's training set constists of. [default: 10000]"

  pvalue_filter_percentile:
    type: int?
    inputBinding:
      prefix: "--par"
      position: 28
    doc: "Percentile for p-value postprocessing filter. [default: 1]"

  poisson:
    type: boolean?
    inputBinding:
      prefix: "--poisson"
      position: 29
    doc: "Use binomial distribution as emmission"

  single_strand:
    type: boolean?
    inputBinding:
      prefix: "--single-strand"
      position: 30
    doc: "Allow single strand BAM file as input"

  m_threshold:
    type: int?
    inputBinding:
      prefix: "--m_threshold"
      position: 31
    doc: "Define the M threshold of percentile for training TMM. [default: 80]"

  a_threshold:
    type: int?
    inputBinding:
      prefix: "--a_threshold"
      position: 32
    doc: "Define the A threshold of percentile for training TMM. [default: 95]"

  remove_duplicates:
    type: boolean?
    inputBinding:
      prefix: "--rmdup"
      position: 33
    doc: "Remove the duplicate reads"


outputs:

  diffpeaks_bed_file:
    type: File?
    outputBinding:
      glob: "*[!-uncor]-diffpeaks.bed"

  diffpeaks_narrowpeak_file:
    type: File?
    outputBinding:
      glob: "*[!-uncor]-diffpeaks.narrowPeak"

  uncor_diffpeaks_bed_file:
    type: File?
    outputBinding:
      glob: "*-uncor-diffpeaks.bed"

  uncor_diffpeaks_narrowpeak_file:
    type: File?
    outputBinding:
      glob: "*-uncor-diffpeaks.narrowPeak"

  cond_1_bigwig_file:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: "*-s1-rep*.bw"

  cond_2_bigwig_file:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: "*-s2-rep*.bw"

  setup_info_file:
    type: File
    outputBinding:
      glob: "*-setup.info"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]
stderr: thor_stderr.log
stdout: thor_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "THOR - differential peak calling of ChIP-seq signals with replicates"
s:name: "THOR - differential peak calling of ChIP-seq signals with replicates"
s:alternateName: "THOR is an HMM-based approach to detect and analyze differential peaks in two sets of ChIP-seq data from distinct biological conditions with replicates"

s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Configuration file is autogenerated based on the bambai_pair_cond_1, bambai_pair_cond_2 and chrom_length_file inputs.
  The following parameters in a configuration file are skipped: genome, inputs1, inputs2.
  The following arguments are skipped: --report (tool fails to execute)
