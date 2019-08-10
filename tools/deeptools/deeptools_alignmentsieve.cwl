cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        if (inputs.output_filename == ""){
            let root = inputs.bambai_pair.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.bambai_pair.basename+".bam":root+".bam";
        } else {
            return inputs.output_filename;
        }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deeptools:v0.0.1


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
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/deeptools_metadata.yaml

s:name: "deeptools_alignmentsieve"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/deeptools/deeptools_alignmentsieve.cwl
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
  For BAM files only. Only selected parameters are implemented.

s:about: |
  usage: Example usage: alignmentSieve.py -b sample1.bam -o sample1.filtered.bam --minMappingQuality 10 --filterMetrics log.txt

  This tool filters alignments in a BAM/CRAM file according the the specified parameters. It can optionally output to BEDPE format.

  optional arguments:
    -h, --help            show this help message and exit

  Required arguments:
    --bam FILE1, -b FILE1
                          An indexed BAM file.
    --outFile OUTFILE, -o OUTFILE
                          The file to write results to. These are the alignments
                          or fragments that pass the filtering criteria.

  General arguments:
    --numberOfProcessors INT, -p INT
                          Number of processors to use. Type "max/2" to use half
                          the maximum number of processors or "max" to use all
                          available processors.
    --filterMetrics FILE.log
                          The number of entries in total and filtered are saved
                          to this file
    --filteredOutReads filtered.bam
                          If desired, all reads NOT passing the filtering
                          criteria can be written to this file.
    --label sample1, -l sample1
                          User defined label instead of the default label (file
                          name).
    --smartLabels         Instead of manually specifying a labels for the input
                          file, this causes deepTools to use the file name after
                          removing the path and extension.
    --verbose, -v         Set to see processing messages.
    --version             show program's version number and exit
    --shift SHIFT [SHIFT ...]
                          Shift the left and right end of a read (for BAM files)
                          or a fragment (for BED files). A positive value shift
                          an end to the right (on the + strand) and a negative
                          value shifts a fragment to the left. Either 2 or 4
                          integers can be provided. For example, "2 -3" will
                          shift the left-most fragment end two bases to the
                          right and the right-most end 3 bases to the left. If 4
                          integers are provided, then the first and last two
                          refer to fragments whose read 1 is on the left or
                          right, respectively. Consequently, it is possible to
                          take strand into consideration for strand-specific
                          protocols. A fragment whose length falls below 1 due
                          to shifting will not be written to the output. See the
                          online documentation for graphical examples. Note that
                          non-properly-paired reads will be filtered.
    --ATACshift           Shift the produced BAM file or BEDPE regions as
                          commonly done for ATAC-seq. This is equivalent to
                          --shift 4 -5 5 -4.

  Output arguments:
    --BED                 Instead of producing BAM files, write output in BEDPE
                          format (as defined by MACS2). Note that only
                          reads/fragments passing filtering criterion are
                          written in BEDPE format.

  Optional arguments:
    --filterRNAstrand {forward,reverse}
                          Selects RNA-seq reads (single-end or paired-end) in
                          the given strand.
    --ignoreDuplicates    If set, reads that have the same orientation and start
                          position will be considered only once. If reads are
                          paired, the mate's position also has to coincide to
                          ignore a read.
    --minMappingQuality INT
                          If set, only reads that have a mapping quality score
                          of at least this are considered.
    --samFlagInclude INT  Include reads based on the SAM flag. For example, to
                          get only reads that are the first mate, use a flag of
                          64. This is useful to count properly paired reads only
                          once, as otherwise the second mate will be also
                          considered for the coverage.
    --samFlagExclude INT  Exclude reads based on the SAM flag. For example, to
                          get only reads that map to the forward strand, use
                          --samFlagExclude 16, where 16 is the SAM flag for
                          reads that map to the reverse strand.
    --blackListFileName BED file [BED file ...], -bl BED file [BED file ...]
                          A BED or GTF file containing regions that should be
                          excluded from all analyses. Currently this works by
                          rejecting genomic chunks that happen to overlap an
                          entry. Consequently, for BAM files, if a read
                          partially overlaps a blacklisted region or a fragment
                          spans over it, then the read/fragment might still be
                          considered. Please note that you should adjust the
                          effective genome size, if relevant.
    --minFragmentLength INT
                          The minimum fragment length needed for read/pair
                          inclusion. This option is primarily useful in ATACseq
                          experiments, for filtering mono- or di-nucleosome
                          fragments.
    --maxFragmentLength INT
                          The maximum fragment length needed for read/pair
                          inclusion.
