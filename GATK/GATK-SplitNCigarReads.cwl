#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

baseCommand:
- gatk
- SplitNCigarReads
doc: |-
  Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).

   Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements).
   The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N
   (including the Ns) is hard clipped and so on for the rest of the new reads. Used for post-processing RNA reads aligned against the full reference.

   <h3>Input</h3>
    <p>
  	    BAM file
    </p>


   <h3>Output</h3>
    <p>
        BAM file with reads split at N CIGAR elements and CIGAR strings updated.
    </p>
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
    - |
      /**
       * File of functions to be added to cwl files
       */

      function generateGATK4BooleanValue(){
          /**
           * Boolean types in GATK 4 are expressed on the command line as --<PREFIX> "true"/"false",
           * so patch here
           */
          if(self === true || self === false){
              return self.toString()
          }

          return self;
      }

      function applyTagsToArgument(prefix, tags){
          /**
           * Function to be used in the field valueFrom of File objects to add gatk tags.
           */

          if(!self){
              return null;
          }
          else if(!tags){
              return generateArrayCmd(prefix);
          }
          else{
              function addTagToArgument(tagObject, argument){
                  var allTags = Array.isArray(tagObject) ? tagObject.join(",") : tagObject;

                  return [prefix + ":" + allTags, argument];
              }

              if(Array.isArray(self)){
                  if(!Array.isArray(tags) || self.length !== tags.length){
                      throw new TypeError("Argument '" + prefix + "' tag field is invalid");
                  }

                  var value = self.map(function(element, i) {
                      return addTagToArgument(tags[i], element);
                  }).reduce(function(a, b){return a.concat(b)})

                  return value;
              }
              else{
                  return addTagToArgument(tags, self);
              }
          }
      }

      function generateArrayCmd(prefix){
          /**
           * Function to be used in the field valueFrom of array objects, so that arrays are optional
           * and prefixes are handled properly.
           *
           * The issue that this solves is documented here:
           * https://www.biostars.org/p/258414/#260140
           */
          if(!self){
              return null;
          }

          if(!Array.isArray(self)){
              self = [self];
          }

          var output = [];
          self.forEach(function(element) {
              output.push(prefix);
              output.push(element);
          })

          return output;
      }

      /* Polyfill String.endsWith (it was introduced in ES6, but CWL 1.0 only supports ES5) */
      String.prototype.endsWith = String.prototype.endsWith || function(suffix) {
          return this.indexOf(suffix, this.length - suffix.length) >= 0;
      };
inputs:
- doc: Reference sequence file [synonymous with -R]
  id: reference
  type: File
  inputBinding:
    prefix: --reference
  secondaryFiles:
  - .fai
  - ^.dict
- doc: BAM/SAM/CRAM file containing reads [synonymous with -I]
  id: reads
  type:
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    prefix: --input
  secondaryFiles: $(self.basename)$(self.nameext.replace('m','i'))?
- doc: Write output to this BAM filename [synonymous with -O]
  id: output_filename
  type: string
  inputBinding:
    prefix: --output
- doc: Read filters to be applied before analysis [synonymous with -RF]
  id: read_filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read-filter"))
- doc: Threshold number of ambiguous bases. If null, uses threshold fraction; otherwise,
    overrides threshold fraction.
  id: ambig-filter-bases
  type: int?
  inputBinding:
    prefix: --ambig-filter-bases
- doc: Threshold fraction of ambiguous bases
  id: ambig-filter-frac
  type: double?
  inputBinding:
    prefix: --ambig-filter-frac
- doc: Maximum length of fragment (insert size)
  id: max-fragment-length
  type: int?
  inputBinding:
    prefix: --max-fragment-length
- doc: Minimum length of fragment (insert size)
  id: min-fragment-length
  type: int?
  inputBinding:
    prefix: --min-fragment-length
- doc: One or more genomic intervals to keep
  id: keep-intervals
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--keep-intervals"))
- doc: Name of the library to keep
  id: library
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--library"))
- doc: Maximum mapping quality to keep (inclusive)
  id: maximum-mapping-quality
  type: int?
  inputBinding:
    prefix: --maximum-mapping-quality
- doc: Minimum mapping quality to keep (inclusive)
  id: minimum-mapping-quality
  type: int?
  inputBinding:
    prefix: --minimum-mapping-quality
- doc: Minimum start location difference at which mapped mates are considered distant
  id: mate-too-distant-length
  type: int?
  inputBinding:
    prefix: --mate-too-distant-length
- doc: Allow a read to be filtered out based on having only 1 soft-clipped block.
    By default, both ends must have a soft-clipped block, setting this flag requires
    only 1 soft-clipped block
  id: dont-require-soft-clips-both-ends
  type: boolean?
  inputBinding:
    prefix: --dont-require-soft-clips-both-ends
    valueFrom: $(generateGATK4BooleanValue())
- doc: Minimum number of aligned bases
  id: filter-too-short
  type: int?
  inputBinding:
    prefix: --filter-too-short
- doc: Platform attribute (PL) to match
  id: platform-filter-name
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--platform-filter-name"))
- doc: Platform unit (PU) to filter out
  id: black-listed-lanes
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--black-listed-lanes"))
- doc: A read group filter expression in the form "attribute:value", where "attribute"
    is a two character read group attribute such as "RG" or "PU".
  id: read-group-black-list
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read-group-black-list"))
- doc: The name of the read group to keep
  id: keep-read-group
  type: string?
  inputBinding:
    prefix: --keep-read-group
- doc: Keep only reads with length at most equal to the specified value
  id: max-read-length
  type: int?
  inputBinding:
    prefix: --max-read-length
- doc: Keep only reads with length at least equal to the specified value
  id: min-read-length
  type: int?
  inputBinding:
    prefix: --min-read-length
- doc: Keep only reads with this read name
  id: read-name
  type: string?
  inputBinding:
    prefix: --read-name
- doc: Keep only reads on the reverse strand
  id: keep-reverse-strand-only
  type: boolean?
  inputBinding:
    prefix: --keep-reverse-strand-only
    valueFrom: $(generateGATK4BooleanValue())
- doc: The name of the sample(s) to keep, filtering out all others
  id: sample
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--sample"))
- doc: Inverts the results from this filter, causing all variants that would pass
    to fail and visa-versa.
  id: invert-soft-clip-ratio-filter
  type: boolean?
  inputBinding:
    prefix: --invert-soft-clip-ratio-filter
    valueFrom: $(generateGATK4BooleanValue())
- doc: Threshold ratio of soft clipped bases (leading / trailing the cigar string)
    to total bases in read for read to be filtered.
  id: soft-clipped-leading-trailing-ratio
  type: double?
  inputBinding:
    prefix: --soft-clipped-leading-trailing-ratio
- doc: Threshold ratio of soft clipped bases (anywhere in the cigar string) to total
    bases in read for read to be filtered.
  id: soft-clipped-ratio-threshold
  type: double?
  inputBinding:
    prefix: --soft-clipped-ratio-threshold
- doc: If true, adds a PG tag to created SAM/BAM/CRAM files.
  id: add-output-sam-program-record
  type: boolean?
  inputBinding:
    prefix: --add-output-sam-program-record
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, adds a command line header line to created VCF files.
  id: add-output-vcf-command-line
  type: boolean?
  inputBinding:
    prefix: --add-output-vcf-command-line
    valueFrom: $(generateGATK4BooleanValue())
- doc: read one or more arguments files and add them to the command line
  id: arguments_file
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--arguments_file", inputs['arguments_file_tags']))
- doc: A argument to set the tags of 'arguments_file'
  id: arguments_file_tags
  type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
- doc: Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer
    if unset. [synonymous with -CIPB]
  id: cloud-index-prefetch-buffer
  type: int?
  inputBinding:
    prefix: --cloud-index-prefetch-buffer
- doc: Size of the cloud-only prefetch buffer (in MB; 0 to disable). [synonymous with
    -CPB]
  id: cloud-prefetch-buffer
  type: int?
  inputBinding:
    prefix: --cloud-prefetch-buffer
- doc: If true, create a BAM/CRAM index when writing a coordinate-sorted BAM/CRAM
    file. [synonymous with -OBI]
  id: create-output-bam-index
  type: boolean?
  inputBinding:
    prefix: --create-output-bam-index
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, create a MD5 digest for any BAM/SAM/CRAM file created [synonymous
    with -OBM]
  id: create-output-bam-md5
  type: boolean?
  inputBinding:
    prefix: --create-output-bam-md5
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, create a VCF index when writing a coordinate-sorted VCF file. [synonymous
    with -OVI]
  id: create-output-variant-index
  type: boolean?
  inputBinding:
    prefix: --create-output-variant-index
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, create a a MD5 digest any VCF file created. [synonymous with -OVM]
  id: create-output-variant-md5
  type: boolean?
  inputBinding:
    prefix: --create-output-variant-md5
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, don't cache bam indexes, this will reduce memory requirements but
    may harm performance if many intervals are specified.  Caching is automatically
    disabled if there are no intervals specified. [synonymous with -DBIC]
  id: disable-bam-index-caching
  type: boolean?
  inputBinding:
    prefix: --disable-bam-index-caching
    valueFrom: $(generateGATK4BooleanValue())
- doc: Read filters to be disabled before analysis [synonymous with -DF]
  id: disable-read-filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--disable-read-filter"))
- doc: If specified, do not check the sequence dictionaries from our inputs for compatibility.
    Use at your own risk!
  id: disable-sequence-dictionary-validation
  type: boolean?
  inputBinding:
    prefix: --disable-sequence-dictionary-validation
    valueFrom: $(generateGATK4BooleanValue())
- doc: 'Disable all tool default read filters (WARNING: many tools will not function
    correctly without their default read filters on)'
  id: disable-tool-default-read-filters
  type: boolean?
  inputBinding:
    prefix: --disable-tool-default-read-filters
    valueFrom: $(generateGATK4BooleanValue())
- doc: do not have the walker soft-clip overhanging sections of the reads
  id: do-not-fix-overhangs
  type: boolean?
  inputBinding:
    prefix: --do-not-fix-overhangs
    valueFrom: $(generateGATK4BooleanValue())
- doc: One or more genomic intervals to exclude from processing [synonymous with -XL]
  id: exclude-intervals
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--exclude-intervals"))
- doc: A configuration file to use with the GATK.
  id: gatk-config-file
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--gatk-config-file", inputs['gatk-config-file_tags']))
- doc: A argument to set the tags of 'gatk-config-file'
  id: gatk-config-file_tags
  type:
  - 'null'
  - string
  - string[]
- doc: If the GCS bucket channel errors out, how many times it will attempt to re-initiate
    the connection [synonymous with -gcs-retries]
  id: gcs-max-retries
  type: int?
  inputBinding:
    prefix: --gcs-max-retries
- doc: Project to bill when accessing "requester pays" buckets. If unset, these buckets
    cannot be accessed.
  id: gcs-project-for-requester-pays
  type: string?
  inputBinding:
    prefix: --gcs-project-for-requester-pays
- doc: A argument to set the tags of 'input'
  id: input_tags
  type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
- doc: Amount of padding (in bp) to add to each interval you are excluding. [synonymous
    with -ixp]
  id: interval-exclusion-padding
  type: int?
  inputBinding:
    prefix: --interval-exclusion-padding
- doc: Interval merging rule for abutting intervals [synonymous with -imr]
  id: interval-merging-rule
  type:
  - 'null'
  - type: enum
    symbols:
    - ALL
    - OVERLAPPING_ONLY
  inputBinding:
    prefix: --interval-merging-rule
- doc: Amount of padding (in bp) to add to each interval you are including. [synonymous
    with -ip]
  id: interval-padding
  type: int?
  inputBinding:
    prefix: --interval-padding
- doc: Set merging approach to use for combining interval inputs [synonymous with
    -isr]
  id: interval-set-rule
  type:
  - 'null'
  - type: enum
    symbols:
    - UNION
    - INTERSECTION
  inputBinding:
    prefix: --interval-set-rule
- doc: One or more genomic intervals over which to operate [synonymous with -L]
  id: intervals
  type:
  - 'null'
  - type: array
    items:
    - File
    - string
    inputBinding:
      valueFrom: $(null)
  - File
  - string
  inputBinding:
    valueFrom: $(applyTagsToArgument("--intervals", inputs['intervals_tags']))
- doc: A argument to set the tags of 'intervals'
  id: intervals_tags
  type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
- doc: Lenient processing of VCF files [synonymous with -LE]
  id: lenient
  type: boolean?
  inputBinding:
    prefix: --lenient
    valueFrom: $(generateGATK4BooleanValue())
- doc: max number of bases allowed in the overhang
  id: max-bases-in-overhang
  type: int?
  inputBinding:
    prefix: --max-bases-in-overhang
- doc: max number of mismatches allowed in the overhang
  id: max-mismatches-in-overhang
  type: int?
  inputBinding:
    prefix: --max-mismatches-in-overhang
- doc: max reads allowed to be kept in memory at a time by the BAM writer
  id: max-reads-in-memory
  type: int?
  inputBinding:
    prefix: --max-reads-in-memory
- doc: have the walker split secondary alignments (will still repair MC tag without
    it)
  id: process-secondary-alignments
  type: boolean?
  inputBinding:
    prefix: --process-secondary-alignments
    valueFrom: $(generateGATK4BooleanValue())
- doc: Whether to suppress job-summary info on System.err.
  id: QUIET
  type: boolean?
  inputBinding:
    prefix: --QUIET
    valueFrom: $(generateGATK4BooleanValue())
- doc: Indices to use for the read inputs. If specified, an index must be provided
    for every read input and in the same order as the read inputs. If this argument
    is not specified, the path to the index for each input will be inferred automatically.
  id: read-index
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read-index"))
- doc: Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The
    default stringency value SILENT can improve performance when processing a BAM
    file in which variable-length data (read, qualities, tags) do not otherwise need
    to be decoded. [synonymous with -VS]
  id: read-validation-stringency
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - LENIENT
    - SILENT
  inputBinding:
    prefix: --read-validation-stringency
- doc: refactor cigar string with NDN elements to one element [synonymous with -fixNDN]
  id: refactor-cigar-string
  type: boolean?
  inputBinding:
    prefix: --refactor-cigar-string
    valueFrom: $(generateGATK4BooleanValue())
- doc: A argument to set the tags of 'reference'
  id: reference_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Output traversal statistics every time this many seconds elapse
  id: seconds-between-progress-updates
  type: double?
  inputBinding:
    prefix: --seconds-between-progress-updates
- doc: Use the given sequence dictionary as the master/canonical sequence dictionary.  Must
    be a .dict file.
  id: sequence-dictionary
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--sequence-dictionary", inputs['sequence-dictionary_tags']))
- doc: A argument to set the tags of 'sequence-dictionary'
  id: sequence-dictionary_tags
  type:
  - 'null'
  - string
  - string[]
- doc: display hidden arguments
  id: showHidden
  type: boolean?
  inputBinding:
    prefix: --showHidden
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, don't emit genotype fields when writing vcf file output.
  id: sites-only-vcf-output
  type: boolean?
  inputBinding:
    prefix: --sites-only-vcf-output
    valueFrom: $(generateGATK4BooleanValue())
- doc: skip the 255 -> 60 MQ read transform [synonymous with -skip-mq-transform]
  id: skip-mapping-quality-transform
  type: boolean?
  inputBinding:
    prefix: --skip-mapping-quality-transform
    valueFrom: $(generateGATK4BooleanValue())
- doc: Temp directory to use.
  id: tmp-dir
  type: string?
  inputBinding:
    prefix: --tmp-dir
- doc: Whether to use the JdkDeflater (as opposed to IntelDeflater) [synonymous with
    -jdk-deflater]
  id: use-jdk-deflater
  type: boolean?
  inputBinding:
    prefix: --use-jdk-deflater
    valueFrom: $(generateGATK4BooleanValue())
- doc: Whether to use the JdkInflater (as opposed to IntelInflater) [synonymous with
    -jdk-inflater]
  id: use-jdk-inflater
  type: boolean?
  inputBinding:
    prefix: --use-jdk-inflater
    valueFrom: $(generateGATK4BooleanValue())
- doc: Control verbosity of logging.
  id: verbosity
  type:
  - 'null'
  - type: enum
    symbols:
    - ERROR
    - WARNING
    - INFO
    - DEBUG
  inputBinding:
    prefix: --verbosity
- doc: display the version number for this tool
  id: version
  type: boolean?
  inputBinding:
    prefix: --version
    valueFrom: $(generateGATK4BooleanValue())
outputs:
- id: output
  doc: Output file from corresponding to the input argument output-filename
  type: File
  outputBinding:
    glob: $(inputs.output_filename)
  secondaryFiles:
  - "$(inputs['create-output-bam-index']? self.basename + self.nameext.replace('m',\
    \ 'i') : [])"
  - "$(inputs['create-output-bam-md5']? self.basename + '.md5' : [])"
