#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        if (inputs.output_filename == ""){
            var root = inputs.bambai_pair.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.bambai_pair.basename+".bam":root+".bam";
        } else {
            return inputs.output_filename;
        }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deeptools:v0.0.1
- class: SoftwareRequirement
  packages:
    deeptools:
      specs: [ "http://identifiers.org/biotools/deeptools" ]
      version: [ "0.0.1" ]

inputs:

  bambai_pair:
    type: File
    secondaryFiles: $(self.basename+".bai")  # due to bug in cwltool==1.0.20190621234233
    inputBinding:
      position: 5
      prefix: "--bam"
    doc: "An indexed BAM file"

  min_fragment_length:
    type: int?
    inputBinding:
      position: 6
      prefix: "--minFragmentLength"
    doc: "The minimum fragment length needed for read/pair inclusion"

  max_fragment_length:
    type: int?
    inputBinding:
      position: 7
      prefix: "--maxFragmentLength"
    doc: "The maximum fragment length needed for read/pair inclusion"

  min_mapping_quality:
    type: int?
    inputBinding:
      position: 8
      prefix: "--minMappingQuality"
    doc: "If set, only reads that have a mapping quality score of at least this are considered"

  ignore_duplicates:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "--ignoreDuplicates"
    doc: |
      If set, reads that have the same orientation and start position will be considered only once.
      If reads are paired, the mate’s position also has to coincide to ignore a read

  shift:
    type:
      - "null"
      - int[]
    inputBinding:
      position: 10
      prefix: "--shift"
    doc: "Shift the left and right end of a read for BAM files"

  atac_shift:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--ATACshift"
    doc: "Shift commonly done for ATAC-seq. This is equivalent to –shift 4 -5 5 -4"

  blacklisted_regions:
    type: File?
    inputBinding:
      position: 12
      prefix: "--blackListFileName"
    doc: "A BED or GTF file containing regions that should be excluded from all analyses"

  output_filename:
    type: string?
    inputBinding:
      position: 13
      prefix: "--outFile"
      valueFrom: $(default_output_filename())
    default: ""
    doc: "The file to write results to. These are the alignments or fragments that pass the filtering criteria"

  threads:
    type: int?
    inputBinding:
      position: 14
      prefix: "--numberOfProcessors"
    doc: "Number of processors to use"

outputs:

  filtered_bam_file:
    type: File
    outputBinding:
      glob: "*.bam"
    doc: "Filtered BAM file"

  alignmentsieve_log:
    type: File
    outputBinding:
      glob: "*.log"
    doc: "alignmentSieve log"

baseCommand: ["alignmentSieve", "--filterMetrics", "alignmentsieve.log"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "deeptools_alignmentsieve"
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
  For BAM files only. Only selected parameters are implemented.