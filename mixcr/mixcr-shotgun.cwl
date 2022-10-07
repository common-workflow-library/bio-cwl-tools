#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1


requirements:
- class: DockerRequirement
  dockerPull: milaboratory/mixcr:3.0.13
- class: ResourceRequirement
  ramMin: 15258
  coresMin: 4
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_prefix = function() {
          if (inputs.output_prefix == ""){
            if ( Array.isArray(inputs.fastq_file) ){
              var root = inputs.fastq_file[0].basename.split('.').slice(0,-1).join('.');
              return (root == "")?inputs.fastq_file[0].basename:root;
            } else {
              var root = inputs.fastq_file.basename.split('.').slice(0,-1).join('.');
              return (root == "")?inputs.fastq_file.basename:root;
            }
          } else {
            return inputs.output_prefix;
          }
        };


inputs:

  species:
    type:
    - type: enum
      symbols:
      - "hsa"
      - "mmu"
      - "rat"
    inputBinding:
      prefix: "--species"
      position: 5
    doc: |
      Species (organism), as specified in library file or taxon id

  starting_material:
    type:
    - type: enum
      symbols:
      - "rna"
      - "dna"
    inputBinding:
      prefix: "--starting-material"
      position: 6
    doc: |
      Type of starting material. This affects which part of reference V
      gene sequences will be used for alignment, with or without intron

  receptor_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "xcr"
      - "tcr"
      - "bcr"
      - "tra"
      - "trb"
      - "trg"
      - "trd"
      - "igh"
      - "igk"
      - "igl"
    inputBinding:
      prefix: "--receptor-type"
      position: 7
    doc: |
      Dedicated receptor type for analysis. By default, all T- and B-cell
      receptor chains are analyzed. MiXCR has special aligner kAligner2,
      which is used when B-cell receptor type is selected. Possible values
      for --receptor-type are: xcr (all chains), tcr, bcr, tra, trb, trg,
      trd, igh, igk, igl.
      Default: xcr (all chains)

  contig_assembly:
    type: boolean?
    inputBinding:
      prefix: "--contig-assembly"
      position: 8
    doc: |
      Whether to assemble full receptor sequences (assembleContigs).
      This option may slow down the computation.

  impute_germline_on_export:
    type: boolean?
    inputBinding:
      prefix: "--impute-germline-on-export"
      position: 9
    doc: |
      Use germline segments (printed with lowercase letters) for
      uncovered gene features.

  only_productive:
    type: boolean?
    inputBinding:
      prefix: "--only-productive"
      position: 10
    doc: |
      Filter out-of-frame sequences and clonotypes with stop-codons
      in clonal sequence export

  assemble_partial_rounds:
    type: int?
    inputBinding:
      prefix: "--assemble-partial-rounds"
      position: 11
    doc: |
      Number of consequent assemblePartial executions (assembles overlapping
      fragmented sequencing reads into long-enough CDR3 containing contigs)
      Default: 2

  do_not_extend_alignments:
    type: boolean?
    inputBinding:
      prefix: "--do-not-extend-alignments"
      position: 12
    doc: |
      Do not perform extension of good TCR alignments

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
      position: 13

  fastq_file:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 14

  output_prefix:
    type: string?
    default: ""
    inputBinding:
      valueFrom: $(default_output_prefix())
      position: 15


outputs:

  alignment_file:
    type: File
    outputBinding:
      glob: $(default_output_prefix() + ".vdjca")  

  assembled_clonotypes_file:
    type: File
    outputBinding:
      glob: $(default_output_prefix() + ".clns")  

  all_clonotypes_file:
    type: File?
    outputBinding:
      glob: $(default_output_prefix() + ".clonotypes.ALL.txt")

  other_clonotypes_files:
    type:
      - "null"
      - type: array
        items: File
    outputBinding:
      glob: $(default_output_prefix() + ".clonotypes.[!A]*.txt")

  report_file:
    type: File
    outputBinding:
      glob: $(default_output_prefix() + ".report")


baseCommand: ["mixcr", "analyze", "shotgun"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "MiXCR Analyze Shotgun"
s:name: "MiXCR Analyze Shotgun"
s:alternateName: "Runs MiXCR pipeline for the analysis of non-enriched RNA-seq and non-targeted genomic data"

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
  MiXCR Analyze Shotgun
  =====================

  MiXCR is a universal framework that processes big immunome data from raw sequences to quantitated clonotypes.
  MiXCR efficiently handles paired- and single-end reads, considers sequence quality, corrects PCR errors and
  identifies germline hypermutations. The software supports both partial- and full-length profiling and employs
  all available RNA or DNA information, including sequences upstream of V and downstream of J gene segments.

  The command analyze shotgun implements the pipeline for the analysis of non-enriched RNA-seq and non-targeted
  genomic data. The pipeline includes alignment of raw sequencing reads using align, assembly of overlapping
  fragmented reads using assemblePartial, imputing good TCR alignments using extend, assembly of aligned sequences
  into clonotypes using assemble and exporting the resulting clonotypes into tab-delimited file using export.
  Optionally, it also assembles full receptor sequences using assembleContigs.

  For more details refer to https://mixcr.readthedocs.io/en/master/rnaseq.html#ref-rna-seq
