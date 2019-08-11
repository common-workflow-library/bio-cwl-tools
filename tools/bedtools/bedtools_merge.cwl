#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.bed_file.basename;
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  bed_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "The input BED file must be sorted by chrom, then start"

  max_distance:
    type: int?
    inputBinding:
      position: 6
      prefix: "-d"
    doc: "Maximum distance between features to be merged" 

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"


outputs:

  merged_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Merged BED file"


baseCommand: ["bedtools", "merge"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools_metadata.yaml

s:name: "bedtools_merge"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bedtools/bedtools_merge.cwl
s:codeRepository: https://github.com/common-workflow-library/bio-cwl-tools
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
  Merges features from BED file. Only selected parameters are implemented.

s:about: |
  Usage:   bedtools merge [OPTIONS] -i <bed/gff/vcf>

  Options: 
    -s	Force strandedness.  That is, only merge features
      that are on the same strand.
      - By default, merging is done without respect to strand.

    -S	Force merge for one specific strand only.
      Follow with + or - to force merge from only
      the forward or reverse strand, respectively.
      - By default, merging is done without respect to strand.

    -d	Maximum distance between features allowed for features
      to be merged.
      - Def. 0. That is, overlapping & book-ended features are merged.
      - (INTEGER)
      - Note: negative values enforce the number of b.p. required for overlap.

    -c	Specify columns from the B file to map onto intervals in A.
      Default: 5.
      Multiple columns can be specified in a comma-delimited list.

    -o	Specify the operation that should be applied to -c.
      Valid operations:
          sum, min, max, absmin, absmax,
          mean, median, mode, antimode
          stdev, sstdev
          collapse (i.e., print a delimited list (duplicates allowed)), 
          distinct (i.e., print a delimited list (NO duplicates allowed)), 
          distinct_sort_num (as distinct, sorted numerically, ascending),
          distinct_sort_num_desc (as distinct, sorted numerically, desscending),
          distinct_only (delimited list of only unique values),
          count
          count_distinct (i.e., a count of the unique values in the column), 
          first (i.e., just the first value in the column), 
          last (i.e., just the last value in the column), 
      Default: sum
      Multiple operations can be specified in a comma-delimited list.

      If there is only column, but multiple operations, all operations will be
      applied on that column. Likewise, if there is only one operation, but
      multiple columns, that operation will be applied to all columns.
      Otherwise, the number of columns must match the the number of operations,
      and will be applied in respective order.
      E.g., "-c 5,4,6 -o sum,mean,count" will give the sum of column 5,
      the mean of column 4, and the count of column 6.
      The order of output columns will match the ordering given in the command.


    -delim	Specify a custom delimiter for the collapse operations.
      - Example: -delim "|"
      - Default: ",".

    -prec	Sets the decimal precision for output (Default: 5)

    -bed	If using BAM input, write output as BED.

    -header	Print the header from the A file prior to results.

    -nobuf	Disable buffered output. Using this option will cause each line
      of output to be printed as it is generated, rather than saved
      in a buffer. This will make printing large output files 
      noticeably slower, but can be useful in conjunction with
      other software tools and scripts that need to process one
      line of bedtools output at a time.

    -iobuf	Specify amount of memory to use for input buffer.
      Takes an integer argument. Optional suffixes K/M/G supported.
      Note: currently has no effect with compressed files.

  Notes: 
    (1) The input file (-i) file must be sorted by chrom, then start.
