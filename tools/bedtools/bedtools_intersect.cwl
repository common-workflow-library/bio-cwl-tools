#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            return inputs.file_a.basename;
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  file_a:
    type: File
    inputBinding:
      position: 5
      prefix: "-a"
    doc: "BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps"

  file_b:
    type: File
    inputBinding:
      position: 6
      prefix: "-b"
    doc: "BAM/BED/GFF/VCF file B. Each feature in A is compared to B in search of overlaps"

  count:
    type: boolean?
    inputBinding:
      position: 7
      prefix: "-c"
    doc: "For each entry in A, report the number of hits in B. Reports 0 for A entries that have no overlap with B" 

  no_overlaps:
    type: boolean?
    inputBinding:
      position: 8
      prefix: "-v"
    doc: "Only report those entries in A that have _no overlaps_ with B" 

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"


outputs:

  intersected_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Intersected BED file"


baseCommand: ["bedtools", "intersect"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools_metadata.yaml

s:name: "bedtools_intersect"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bedtools/bedtools_intersect.cwl
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
  Intersect features from A and B file. Only selected parameters are implemented.

s:about: |
  Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>

    Note: -b may be followed with multiple databases and/or 
    wildcard (*) character(s). 
  Options: 
    -wa	Write the original entry in A for each overlap.

    -wb	Write the original entry in B for each overlap.
      - Useful for knowing _what_ A overlaps. Restricted by -f and -r.

    -loj	Perform a "left outer join". That is, for each feature in A
      report each overlap with B.  If no overlaps are found, 
      report a NULL feature for B.

    -wo	Write the original A and B entries plus the number of base
      pairs of overlap between the two features.
      - Overlaps restricted by -f and -r.
        Only A features with overlap are reported.

    -wao	Write the original A and B entries plus the number of base
      pairs of overlap between the two features.
      - Overlapping features restricted by -f and -r.
        However, A features w/o overlap are also reported
        with a NULL B feature and overlap = 0.

    -u	Write the original A entry _once_ if _any_ overlaps found in B.
      - In other words, just report the fact >=1 hit was found.
      - Overlaps restricted by -f and -r.

    -c	For each entry in A, report the number of overlaps with B.
      - Reports 0 for A entries that have no overlap with B.
      - Overlaps restricted by -f and -r.

    -v	Only report those entries in A that have _no overlaps_ with B.
      - Similar to "grep -v" (an homage).

    -ubam	Write uncompressed BAM output. Default writes compressed BAM.

    -s	Require same strandedness.  That is, only report hits in B
      that overlap A on the _same_ strand.
      - By default, overlaps are reported without respect to strand.

    -S	Require different strandedness.  That is, only report hits in B
      that overlap A on the _opposite_ strand.
      - By default, overlaps are reported without respect to strand.

    -f	Minimum overlap required as a fraction of A.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -F	Minimum overlap required as a fraction of B.
      - Default is 1E-9 (i.e., 1bp).
      - FLOAT (e.g. 0.50)

    -r	Require that the fraction overlap be reciprocal for A AND B.
      - In other words, if -f is 0.90 and -r is used, this requires
        that B overlap 90% of A and A _also_ overlaps 90% of B.

    -e	Require that the minimum fraction be satisfied for A OR B.
      - In other words, if -e is used with -f 0.90 and -F 0.10 this requires
        that either 90% of A is covered OR 10% of  B is covered.
        Without -e, both fractions would have to be satisfied.

    -split	Treat "split" BAM or BED12 entries as distinct BED intervals.

    -g	Provide a genome file to enforce consistent chromosome sort order
      across input files. Only applies when used with -sorted option.

    -nonamecheck	For sorted data, don't throw an error if the file has different naming conventions
        for the same chromosome. ex. "chr1" vs "chr01".

    -sorted	Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input.

    -names	When using multiple databases, provide an alias for each that
      will appear instead of a fileId when also printing the DB record.

    -filenames	When using multiple databases, show each complete filename
        instead of a fileId when also printing the DB record.

    -sortout	When using multiple databases, sort the output DB hits
        for each record.

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
    (1) When a BAM file is used for the A file, the alignment is retained if overlaps exist,
    and exlcuded if an overlap cannot be found.  If multiple overlaps exist, they are not
    reported, as we are only testing for one or more overlaps.
