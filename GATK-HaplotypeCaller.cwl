#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand:
- gatk
- HaplotypeCaller
doc: |-
  Call germline SNPs and indels via local re-assembly of haplotypes

   <p>The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.</p>

   <p>In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).</p>

   <p>In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. For that purpose, use Mutect2 instead.</p>

   <p>Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge for most variant callers,
   on the condition that the input read data has previously been processed according to our recommendations as documented <a href='https://software.broadinstitute.org/gatk/documentation/article?id=4067'>here</a>.</p>

   <h3>How HaplotypeCaller works</h3>

   <br />
   <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4147'>1. Define active regions </a></h4>

   <p>The program determines which regions of the genome it needs to operate on (active regions), based on the presence of
   evidence for variation.

   <br />
   <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4146'>2. Determine haplotypes by assembly of the active region </a></h4>

   <p>For each active region, the program builds a De Bruijn-like graph to reassemble the active region and identifies
   what are the possible haplotypes present in the data. The program then realigns each haplotype against the reference
   haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites. </p>

   <br />
   <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4441'>3. Determine likelihoods of the haplotypes given the read data </a></h4>

   <p>For each active region, the program performs a pairwise alignment of each read against each haplotype using the
   PairHMM algorithm. This produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are
   then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.   </p>

   <br />
   <h4><a href='https://software.broadinstitute.org/gatk/documentation/article?id=4442'>4. Assign sample genotypes </a></h4>

   <p>For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the
   read data to calculate the likelihoods of each genotype per sample given the read data observed for that
   sample. The most likely genotype is then assigned to the sample.    </p>

   <h3>Input</h3>
   <p>
   Input bam file(s) from which to make variant calls
   </p>

   <h3>Output</h3>
   <p>
   Either a VCF or GVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant
   recalibration (Best Practice) or hard-filtering before use in downstream analyses. If using the GVCF workflow, the
   output is a GVCF file that must first be run through GenotypeGVCFs and then filtering before further analysis.
   </p>

   <h3>Caveats</h3>
   <ul>
   <li>We have not yet fully tested the interaction between the GVCF-based calling or the multisample calling and the
   RNAseq-specific functionalities. Use those in combination at your own risk.</li>
   </ul>

   <h3>Special note on ploidy</h3>
   <p>This tool is able to handle many non-diploid use cases; the desired ploidy can be specified using the -ploidy
   argument. Note however that very high ploidies (such as are encountered in large pooled experiments) may cause
   performance challenges including excessive slowness. We are working on resolving these limitations.</p>

   <h3>Additional Notes</h3>
   <ul>
       <li>When working with PCR-free data, be sure to set `-pcr_indel_model NONE` (see argument below).</li>
       <li>When running in `-ERC GVCF` or `-ERC BP_RESOLUTION` modes, the confidence threshold
       is automatically set to 0. This cannot be overridden by the command line. The threshold can be set manually
       to the desired level in the next step of the workflow (GenotypeGVCFs)</li>
       <li>We recommend using a list of intervals to speed up analysis. See <a href='https://software.broadinstitute.org/gatk/documentation/article?id=4133'>this document</a> for details.</li>
   </ul>
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
  SchemaDefRequirement:
    types:
    - type: enum
      name: annotation_type
      symbols:
      - AS_BaseQualityRankSumTest
      - AS_FisherStrand
      - AS_InbreedingCoeff
      - AS_MappingQualityRankSumTest
      - AS_QualByDepth
      - AS_RMSMappingQuality
      - AS_ReadPosRankSumTest
      - AS_StrandOddsRatio
      - AlleleFraction
      - BaseQuality
      - BaseQualityRankSumTest
      - ChromosomeCounts
      - ClippingRankSumTest
      - CountNs
      - Coverage
      - DepthPerAlleleBySample
      - DepthPerSampleHC
      - ExcessHet
      - FisherStrand
      - FragmentLength
      - GenotypeSummaries
      - InbreedingCoeff
      - LikelihoodRankSumTest
      - MappingQuality
      - MappingQualityRankSumTest
      - MappingQualityZero
      - OrientationBiasReadCounts
      - OriginalAlignment
      - PossibleDeNovo
      - QualByDepth
      - RMSMappingQuality
      - ReadPosRankSumTest
      - ReadPosition
      - ReferenceBases
      - SampleList
      - StrandBiasBySample
      - StrandOddsRatio
      - TandemRepeat
      - UniqueAltReadCount
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0
inputs:
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
- doc: Minimum probability for a locus to be considered active.
  id: active-probability-threshold
  type: double?
  inputBinding:
    prefix: --active-probability-threshold
- doc: Use Mutect2's adaptive graph pruning algorithm
  id: adaptive-pruning
  type: boolean?
  inputBinding:
    prefix: --adaptive-pruning
    valueFrom: $(generateGATK4BooleanValue())
- doc: Initial base error rate estimate for adaptive pruning
  id: adaptive-pruning-initial-error-rate
  type: double?
  inputBinding:
    prefix: --adaptive-pruning-initial-error-rate
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
- doc: Annotate all sites with PLs
  id: all-site-pls
  type: boolean?
  inputBinding:
    prefix: --all-site-pls
    valueFrom: $(generateGATK4BooleanValue())
- doc: Likelihood and read-based annotations will only take into consideration reads
    that overlap the variant or any base no further than this distance expressed in
    base pairs
  id: allele-informative-reads-overlap-margin
  type: int?
  inputBinding:
    prefix: --allele-informative-reads-overlap-margin
- doc: The set of alleles to force-call regardless of evidence
  id: alleles
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--alleles", inputs['alleles_tags']))
- doc: A argument to set the tags of 'alleles'
  id: alleles_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Allow graphs that have non-unique kmers in the reference
  id: allow-non-unique-kmers-in-ref
  type: boolean?
  inputBinding:
    prefix: --allow-non-unique-kmers-in-ref
    valueFrom: $(generateGATK4BooleanValue())
- doc: If provided, we will annotate records with the number of alternate alleles
    that were discovered (but not necessarily genotyped) at a given site
  id: annotate-with-num-discovered-alleles
  type: boolean?
  inputBinding:
    prefix: --annotate-with-num-discovered-alleles
    valueFrom: $(generateGATK4BooleanValue())
- doc: One or more specific annotations to add to variant calls [synonymous with -A]
  id: annotation
  type:
  - 'null'
  - annotation_type
  - annotation_type[]
  inputBinding:
    valueFrom: $(generateArrayCmd("--annotation"))
- doc: One or more groups of annotations to apply to variant calls [synonymous with
    -G]
  id: annotation-group
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--annotation-group"))
- doc: One or more specific annotations to exclude from variant calls [synonymous
    with -AX]
  id: annotations-to-exclude
  type:
  - 'null'
  - annotation_type
  - annotation_type[]
  inputBinding:
    valueFrom: $(generateArrayCmd("--annotations-to-exclude"))
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
- doc: Output the assembly region to this IGV formatted file
  id: assembly-region-out-filename
  type: string?
  inputBinding:
    prefix: --assembly-region-out
- doc: Number of additional bases of context to include around each assembly region
  id: assembly-region-padding
  type: int?
  inputBinding:
    prefix: --assembly-region-padding
- doc: File to which assembled haplotypes should be written [synonymous with -bamout]
  id: bam-output-filename
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--bam-output", inputs['bam-output_tags']))
- doc: A argument to set the tags of 'bam-output'
  id: bam-output_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Which haplotypes should be written to the BAM
  id: bam-writer-type
  type:
  - 'null'
  - type: enum
    symbols:
    - ALL_POSSIBLE_HAPLOTYPES
    - CALLED_HAPLOTYPES
  inputBinding:
    prefix: --bam-writer-type
- doc: Base qualities below this threshold will be reduced to the minimum (6)
  id: base-quality-score-threshold
  type: int?
  inputBinding:
    prefix: --base-quality-score-threshold
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
- doc: Comparison VCF file(s) [synonymous with -comp]
  id: comparison
  type:
  - 'null'
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--comparison", inputs['comparison_tags']))
- doc: A argument to set the tags of 'comparison'
  id: comparison_tags
  type:
  - 'null'
  - type: array
    items:
    - string
    - type: array
      items: string
- doc: Tab-separated File containing fraction of contamination in sequencing data
    (per sample) to aggressively remove. Format should be "<SampleID><TAB><Contamination>"
    (Contamination is double) per line; No header. [synonymous with -contamination-file]
  id: contamination-fraction-per-sample-file
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--contamination-fraction-per-sample-file", inputs['contamination-fraction-per-sample-file_tags']))
- doc: A argument to set the tags of 'contamination-fraction-per-sample-file'
  id: contamination-fraction-per-sample-file_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Fraction of contamination in sequencing data (for all samples) to aggressively
    remove [synonymous with -contamination]
  id: contamination-fraction-to-filter
  type: double?
  inputBinding:
    prefix: --contamination-fraction-to-filter
- doc: Undocumented option
  id: correct-overlapping-quality
  type: boolean?
  inputBinding:
    prefix: --correct-overlapping-quality
    valueFrom: $(generateGATK4BooleanValue())
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
- doc: dbSNP file [synonymous with -D]
  id: dbsnp
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--dbsnp", inputs['dbsnp_tags']))
- doc: A argument to set the tags of 'dbsnp'
  id: dbsnp_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Print out verbose debug information about each assembly region [synonymous
    with -debug]
  id: debug-assembly
  type: boolean?
  inputBinding:
    prefix: --debug-assembly
    valueFrom: $(generateGATK4BooleanValue())
- doc: If true, don't cache bam indexes, this will reduce memory requirements but
    may harm performance if many intervals are specified.  Caching is automatically
    disabled if there are no intervals specified. [synonymous with -DBIC]
  id: disable-bam-index-caching
  type: boolean?
  inputBinding:
    prefix: --disable-bam-index-caching
    valueFrom: $(generateGATK4BooleanValue())
- doc: Don't skip calculations in ActiveRegions with no variants
  id: disable-optimizations
  type: boolean?
  inputBinding:
    prefix: --disable-optimizations
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
- doc: Disable all tool default annotations
  id: disable-tool-default-annotations
  type: boolean?
  inputBinding:
    prefix: --disable-tool-default-annotations
    valueFrom: $(generateGATK4BooleanValue())
- doc: 'Disable all tool default read filters (WARNING: many tools will not function
    correctly without their default read filters on)'
  id: disable-tool-default-read-filters
  type: boolean?
  inputBinding:
    prefix: --disable-tool-default-read-filters
    valueFrom: $(generateGATK4BooleanValue())
- doc: Disable physical phasing
  id: do-not-run-physical-phasing
  type: boolean?
  inputBinding:
    prefix: --do-not-run-physical-phasing
    valueFrom: $(generateGATK4BooleanValue())
- doc: Disable iterating over kmer sizes when graph cycles are detected
  id: dont-increase-kmer-sizes-for-cycles
  type: boolean?
  inputBinding:
    prefix: --dont-increase-kmer-sizes-for-cycles
    valueFrom: $(generateGATK4BooleanValue())
- doc: Do not analyze soft clipped bases in the reads
  id: dont-use-soft-clipped-bases
  type: boolean?
  inputBinding:
    prefix: --dont-use-soft-clipped-bases
    valueFrom: $(generateGATK4BooleanValue())
- doc: Mode for emitting reference confidence scores (For Mutect2, this is a BETA
    feature) [synonymous with -ERC]
  id: emit-ref-confidence
  type:
  - 'null'
  - type: enum
    symbols:
    - NONE
    - BP_RESOLUTION
    - GVCF
  inputBinding:
    prefix: --emit-ref-confidence
- doc: Use all possible annotations (not for the faint of heart)
  id: enable-all-annotations
  type: boolean?
  inputBinding:
    prefix: --enable-all-annotations
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
- doc: Output the band lower bound for each GQ block regardless of the data it represents
  id: floor-blocks
  type: boolean?
  inputBinding:
    prefix: --floor-blocks
    valueFrom: $(generateGATK4BooleanValue())
- doc: If provided, all regions will be marked as active
  id: force-active
  type: boolean?
  inputBinding:
    prefix: --force-active
    valueFrom: $(generateGATK4BooleanValue())
- doc: Force-call filtered alleles included in the resource specified by --alleles
    [synonymous with -genotype-filtered-alleles]
  id: force-call-filtered-alleles
  type: boolean?
  inputBinding:
    prefix: --force-call-filtered-alleles
    valueFrom: $(generateGATK4BooleanValue())
- doc: Samples representing the population "founders"
  id: founder-id
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--founder-id"))
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
- doc: Write debug assembly graph information to this file [synonymous with -graph]
  id: graph-output-filename
  type: string?
  inputBinding:
    prefix: --graph-output
- doc: Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100]
    and specified in increasing order) [synonymous with -GQB]
  id: gvcf-gq-bands
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--gvcf-gq-bands"))
- doc: Heterozygosity value used to compute prior likelihoods for any locus.  See
    the GATKDocs for full details on the meaning of this population genetics concept
  id: heterozygosity
  type: double?
  inputBinding:
    prefix: --heterozygosity
- doc: Standard deviation of heterozygosity for SNP and indel calling.
  id: heterozygosity-stdev
  type: double?
  inputBinding:
    prefix: --heterozygosity-stdev
- doc: Heterozygosity for indel calling.  See the GATKDocs for heterozygosity for
    full details on the meaning of this population genetics concept
  id: indel-heterozygosity
  type: double?
  inputBinding:
    prefix: --indel-heterozygosity
- doc: The size of an indel to check for in the reference model
  id: indel-size-to-eliminate-in-ref-model
  type: int?
  inputBinding:
    prefix: --indel-size-to-eliminate-in-ref-model
- doc: BAM/SAM/CRAM file containing reads [synonymous with -I]
  id: input
  type:
  - type: array
    items: File
    inputBinding:
      valueFrom: $(null)
  - File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--input", inputs['input_tags']))
  secondaryFiles: $(self.basename + self.nameext.replace('m','i'))
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
- doc: Kmer size to use in the read threading assembler
  id: kmer-size
  type:
  - 'null'
  - type: array
    items: int
    inputBinding:
      valueFrom: $(null)
  - int
  inputBinding:
    valueFrom: $(generateArrayCmd("--kmer-size"))
- doc: Lenient processing of VCF files [synonymous with -LE]
  id: lenient
  type: boolean?
  inputBinding:
    prefix: --lenient
    valueFrom: $(generateGATK4BooleanValue())
- doc: Maximum number of alternate alleles to genotype
  id: max-alternate-alleles
  type: int?
  inputBinding:
    prefix: --max-alternate-alleles
- doc: Maximum size of an assembly region
  id: max-assembly-region-size
  type: int?
  inputBinding:
    prefix: --max-assembly-region-size
- doc: Maximum number of genotypes to consider at any site
  id: max-genotype-count
  type: int?
  inputBinding:
    prefix: --max-genotype-count
- doc: Two or more phased substitutions separated by this distance or less are merged
    into MNPs. [synonymous with -mnp-dist]
  id: max-mnp-distance
  type: int?
  inputBinding:
    prefix: --max-mnp-distance
- doc: Maximum number of haplotypes to consider for your population
  id: max-num-haplotypes-in-population
  type: int?
  inputBinding:
    prefix: --max-num-haplotypes-in-population
- doc: Upper limit on how many bases away probability mass can be moved around when
    calculating the boundaries between active and inactive assembly regions
  id: max-prob-propagation-distance
  type: int?
  inputBinding:
    prefix: --max-prob-propagation-distance
- doc: Maximum number of reads to retain per alignment start position. Reads above
    this threshold will be downsampled. Set to 0 to disable.
  id: max-reads-per-alignment-start
  type: int?
  inputBinding:
    prefix: --max-reads-per-alignment-start
- doc: Maximum number of variants in graph the adaptive pruner will allow
  id: max-unpruned-variants
  type: int?
  inputBinding:
    prefix: --max-unpruned-variants
- doc: Minimum size of an assembly region
  id: min-assembly-region-size
  type: int?
  inputBinding:
    prefix: --min-assembly-region-size
- doc: Minimum base quality required to consider a base for calling [synonymous with
    -mbq]
  id: min-base-quality-score
  type: int?
  inputBinding:
    prefix: --min-base-quality-score
- doc: Minimum length of a dangling branch to attempt recovery
  id: min-dangling-branch-length
  type: int?
  inputBinding:
    prefix: --min-dangling-branch-length
- doc: Minimum support to not prune paths in the graph
  id: min-pruning
  type: int?
  inputBinding:
    prefix: --min-pruning
- doc: How many threads should a native pairHMM implementation use
  id: native-pair-hmm-threads
  type: int?
  inputBinding:
    prefix: --native-pair-hmm-threads
- doc: use double precision in the native pairHmm. This is slower but matches the
    java implementation better
  id: native-pair-hmm-use-double-precision
  type: boolean?
  inputBinding:
    prefix: --native-pair-hmm-use-double-precision
    valueFrom: $(generateGATK4BooleanValue())
- doc: Number of samples that must pass the minPruning threshold
  id: num-pruning-samples
  type: int?
  inputBinding:
    prefix: --num-pruning-samples
- doc: Number of hom-ref genotypes to infer at sites not present in a panel
  id: num-reference-samples-if-no-call
  type: int?
  inputBinding:
    prefix: --num-reference-samples-if-no-call
- doc: File to which variants should be written [synonymous with -O]
  id: output_filename
  type: string
  inputBinding:
    prefix: --output
- doc: Specifies which type of calls we should output
  id: output_mode
  type:
  - 'null'
  - type: enum
    symbols:
    - EMIT_VARIANTS_ONLY
    - EMIT_ALL_CONFIDENT_SITES
    - EMIT_ALL_ACTIVE_SITES
  inputBinding:
    prefix: --output-mode
- doc: Flat gap continuation penalty for use in the Pair HMM
  id: pair-hmm-gap-continuation-penalty
  type: int?
  inputBinding:
    prefix: --pair-hmm-gap-continuation-penalty
- doc: The PairHMM implementation to use for genotype likelihood calculations [synonymous
    with -pairHMM]
  id: pair-hmm-implementation
  type:
  - 'null'
  - type: enum
    symbols:
    - EXACT
    - ORIGINAL
    - LOGLESS_CACHING
    - AVX_LOGLESS_CACHING
    - AVX_LOGLESS_CACHING_OMP
    - EXPERIMENTAL_FPGA_LOGLESS_CACHING
    - FASTEST_AVAILABLE
  inputBinding:
    prefix: --pair-hmm-implementation
- doc: The PCR indel model to use
  id: pcr-indel-model
  type:
  - 'null'
  - type: enum
    symbols:
    - NONE
    - HOSTILE
    - AGGRESSIVE
    - CONSERVATIVE
  inputBinding:
    prefix: --pcr-indel-model
- doc: Pedigree file for determining the population "founders" [synonymous with -ped]
  id: pedigree
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--pedigree", inputs['pedigree_tags']))
- doc: A argument to set the tags of 'pedigree'
  id: pedigree_tags
  type:
  - 'null'
  - string
  - string[]
- doc: The global assumed mismapping rate for reads
  id: phred-scaled-global-read-mismapping-rate
  type: int?
  inputBinding:
    prefix: --phred-scaled-global-read-mismapping-rate
- doc: Callset to use in calculating genotype priors [synonymous with -population]
  id: population-callset
  type: File?
  inputBinding:
    valueFrom: $(applyTagsToArgument("--population-callset", inputs['population-callset_tags']))
- doc: A argument to set the tags of 'population-callset'
  id: population-callset_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Ln likelihood ratio threshold for adaptive pruning algorithm
  id: pruning-lod-threshold
  type: double?
  inputBinding:
    prefix: --pruning-lod-threshold
- doc: Whether to suppress job-summary info on System.err.
  id: QUIET
  type: boolean?
  inputBinding:
    prefix: --QUIET
    valueFrom: $(generateGATK4BooleanValue())
- doc: Read filters to be applied before analysis [synonymous with -RF]
  id: read-filter
  type:
  - 'null'
  - type: array
    items: string
    inputBinding:
      valueFrom: $(null)
  - string
  inputBinding:
    valueFrom: $(generateArrayCmd("--read-filter"))
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
- doc: Recover all dangling branches
  id: recover-all-dangling-branches
  type: boolean?
  inputBinding:
    prefix: --recover-all-dangling-branches
    valueFrom: $(generateGATK4BooleanValue())
- doc: This argument is deprecated since version 3.3
  id: recover-dangling-heads
  type: boolean?
  inputBinding:
    prefix: --recover-dangling-heads
    valueFrom: $(generateGATK4BooleanValue())
- doc: Reference sequence file [synonymous with -R]
  id: reference
  type: File
  inputBinding:
    valueFrom: $(applyTagsToArgument("--reference", inputs['reference_tags']))
  secondaryFiles:
  - .fai
  - ^.dict
- doc: A argument to set the tags of 'reference'
  id: reference_tags
  type:
  - 'null'
  - string
  - string[]
- doc: Name of single sample to use from a multi-sample bam [synonymous with -ALIAS]
  id: sample-name
  type: string?
  inputBinding:
    prefix: --sample-name
- doc: Ploidy (number of chromosomes) per sample. For pooled data, set to (Number
    of samples in each pool * Sample Ploidy). [synonymous with -ploidy]
  id: sample-ploidy
  type: int?
  inputBinding:
    prefix: --sample-ploidy
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
- doc: Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is
    the right choice
  id: smith-waterman
  type:
  - 'null'
  - type: enum
    symbols:
    - FASTEST_AVAILABLE
    - AVX_ENABLED
    - JAVA
  inputBinding:
    prefix: --smith-waterman
- doc: The minimum phred-scaled confidence threshold at which variants should be called
    [synonymous with -stand-call-conf]
  id: standard-min-confidence-threshold-for-calling
  type: double?
  inputBinding:
    prefix: --standard-min-confidence-threshold-for-calling
- doc: Temp directory to use.
  id: tmp-dir
  type: string?
  inputBinding:
    prefix: --tmp-dir
- doc: Use the contamination-filtered read maps for the purposes of annotating variants
  id: use-filtered-reads-for-annotations
  type: boolean?
  inputBinding:
    prefix: --use-filtered-reads-for-annotations
    valueFrom: $(generateGATK4BooleanValue())
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
- id: assembly-region-out
  doc: Output file from corresponding to the input argument assembly-region-out-filename
  type: File?
  outputBinding:
    glob: $(inputs['assembly-region-out-filename'])
- id: bam-output
  doc: Output file from corresponding to the input argument bam-output-filename
  type: File?
  outputBinding:
    glob: $(inputs['bam-output-filename'])
  secondaryFiles:
  - "$(inputs['create-output-bam-index']? self.basename + self.nameext.replace('m',\
    \ 'i') : [])"
  - "$(inputs['create-output-bam-md5']? self.basename + '.md5' : [])"
- id: graph-output
  doc: Output file from corresponding to the input argument graph-output-filename
  type: File?
  outputBinding:
    glob: $(inputs['graph-output-filename'])
- id: output
  doc: Output file from corresponding to the input argument output-filename
  type: File
  outputBinding:
    glob: $(inputs.output_filename)
  secondaryFiles:
  - "$(inputs['create-output-variant-index']? self.basename + (inputs.output_filename.endsWith('.gz')?\
    \ '.tbi':'.idx') : [])"
  - "$(inputs['create-output-variant-md5']? self.basename + '.md5' : [])"
