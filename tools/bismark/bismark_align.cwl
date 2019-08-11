#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2


inputs:

  indices_folder:
    type: Directory
    inputBinding:
      position: 3
    label: "Bismark indices folder"
    doc: "Path to Bismark generated indices folder"

  fastq_file:
    type: File
    inputBinding:
      position: 4
    label: "FASTQ file"
    doc: "Uncompressed or gzipped FASTQ file, single-end"

  processes:
    type: int?
    inputBinding:
      position: 1
      prefix: "--multicore"
    label: "Number of Bismark instances to run"
    doc: "Set the number of parallel Bismark instances to run concurrently. Each Bismark instance runs four Bowtie2 aligners"

  threads:
    type: int?
    inputBinding:
      position: 2
      prefix: "-p"
    label: "Number of Bowtie2 threads to use"
    doc: "Set the number of threads for each Bowtie2 aligner"


outputs:

  bam_file:
    type: File
    label: "BAM alignment file"
    doc: "Bismark generated BAM alignment file"
    outputBinding:
      glob: "*.bam"

  alignment_report:
    type: File
    label: "Bismark alignment and methylation report"
    doc: "Bismark generated alignment and methylation summary report"
    outputBinding:
      glob: "*.txt"


baseCommand: ["bismark", "--non_directional"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bismark_metadata.yaml

s:name: "bismark_align"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bismark/bismark_align.cwl
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
  Default aligner - Bowtie2.
  Only Single-End supported.
  Parameters used:
  --non_directional
    The sequencing library was constructed in a non strand-specific manner, alignments to all four bisulfite strands
    will be reported. (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary
    to the original strands are merely theoretical and should not exist in reality. Specifying directional alignments
    (which is the default) will only run 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel
    and report these alignments. This is the recommended option for strand-specific libraries).

s:about: |
  DESCRIPTION

  The following is a brief description of command line options and arguments to control the Bismark
  bisulfite mapper and methylation caller. Bismark takes in FastA or FastQ files and aligns the
  reads to a specified bisulfite genome. Sequence reads are transformed into a bisulfite converted forward strand
  version (C->T conversion) or into a bisulfite treated reverse strand (G->A conversion of the forward strand).
  Each of these reads are then aligned to bisulfite treated forward strand index of a reference genome
  (C->T converted) and a bisulfite treated reverse strand index of the genome (G->A conversion of the
  forward strand, by doing this alignments will produce the same positions). These 4 instances of Bowtie (1 or 2)
  are run in parallel. The sequence file(s) are then read in again sequence by sequence to pull out the original
  sequence from the genome and determine if there were any protected C's present or not.

  As of version 0.7.0 Bismark will only run 2 alignment threads for OT and OB in parallel, the 4 strand mode can be
  re-enabled by using --non_directional.

  The final output of Bismark is in SAM format by default. For Bowtie 1 one can alos choose to report the old
  'vanilla' output format, which is a single tab delimited file with all sequences that have a unique best
  alignment to any of the 4 possible strands of a bisulfite PCR product. Both formats are described in more detail below.


  USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}


  ARGUMENTS:

  <genome_folder>          The path to the folder containing the unmodified reference genome
                          as well as the subfolders created by the Bismark_Genome_Preparation
                          script (/Bisulfite_Genome/CT_conversion/ and /Bisulfite_Genome/GA_conversion/).
                          Bismark expects one or more fastA files in this folder (file extension: .fa
                          or .fasta). The path can be relative or absolute. The path may also be set as
                          '--genome_folder /path/to/genome/folder/'.

  -1 <mates1>              Comma-separated list of files containing the #1 mates (filename usually includes
                          "_1"), e.g. flyA_1.fq,flyB_1.fq). Sequences specified with this option must
                          correspond file-for-file and read-for-read with those specified in <mates2>.
                          Reads may be a mix of different lengths. Bismark will produce one mapping result
                          and one report file per paired-end input file pair.

  -2 <mates2>              Comma-separated list of files containing the #2 mates (filename usually includes
                          "_2"), e.g. flyA_1.fq,flyB_1.fq). Sequences specified with this option must
                          correspond file-for-file and read-for-read with those specified in <mates1>.
                          Reads may be a mix of different lengths.

  <singles>                A comma- or space-separated list of files containing the reads to be aligned (e.g.
                          lane1.fq,lane2.fq lane3.fq). Reads may be a mix of different lengths. Bismark will
                          produce one mapping result and one report file per input file.


  OPTIONS:


  Input:

  --se/--single_end <list> Sets single-end mapping mode explicitly giving a list of file names as <list>.
                          The filenames may be provided as a comma [,] or colon [:] separated list.

  -q/--fastq               The query input files (specified as <mate1>,<mate2> or <singles> are FASTQ
                          files (usually having extension .fg or .fastq). This is the default. See also
                          --solexa-quals.

  -f/--fasta               The query input files (specified as <mate1>,<mate2> or <singles> are FASTA
                          files (usually having extensions .fa, .mfa, .fna or similar). All quality values
                          are assumed to be 40 on the Phred scale. FASTA files are expected to contain both
                          the read name and the sequence on a single line (and not spread over several lines).

  -s/--skip <int>          Skip (i.e. do not align) the first <int> reads or read pairs from the input.

  -u/--upto <int>          Only aligns the first <int> reads or read pairs from the input. Default: no limit.

  --phred33-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 33. Default: on.

  --phred64-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 64. Default: off.

  --solexa-quals           Convert FASTQ qualities from solexa-scaled (which can be negative) to phred-scaled
                          (which can't). The formula for conversion is:
                          phred-qual = 10 * log(1 + 10 ** (solexa-qual/10.0)) / log(10). Used with -q. This
                          is usually the right option for use with (unconverted) reads emitted by the GA
                          Pipeline versions prior to 1.3. Works only for Bowtie 1. Default: off.

  --solexa1.3-quals        Same as --phred64-quals. This is usually the right option for use with (unconverted)
                          reads emitted by GA Pipeline version 1.3 or later. Default: off.

  --path_to_bowtie         The full path </../../> to the Bowtie (1 or 2) installation on your system. If not
                          specified it is assumed that Bowtie (1 or 2) is in the PATH.


  Alignment:

  -n/--seedmms <int>       The maximum number of mismatches permitted in the "seed", i.e. the first L base pairs
                          of the read (where L is set with -l/--seedlen). This may be 0, 1, 2 or 3 and the
                          default is 1. This option is only available for Bowtie 1 (for Bowtie 2 see -N).

  -l/--seedlen             The "seed length"; i.e., the number of bases of the high quality end of the read to
                          which the -n ceiling applies. The default is 28. Bowtie (and thus Bismark) is faster for
                          larger values of -l. This option is only available for Bowtie 1 (for Bowtie 2 see -L).

  -e/--maqerr <int>        Maximum permitted total of quality values at all mismatched read positions throughout
                          the entire alignment, not just in the "seed". The default is 70. Like Maq, bowtie rounds
                          quality values to the nearest 10 and saturates at 30. This value is not relevant for
                          Bowtie 2.

  --chunkmbs <int>         The number of megabytes of memory a given thread is given to store path descriptors in
                          --best mode. Best-first search must keep track of many paths at once to ensure it is
                          always extending the path with the lowest cumulative cost. Bowtie tries to minimize the
                          memory impact of the descriptors, but they can still grow very large in some cases. If
                          you receive an error message saying that chunk memory has been exhausted in --best mode,
                          try adjusting this parameter up to dedicate more memory to the descriptors. This value
                          is not relevant for Bowtie 2. Default: 512.

  -I/--minins <int>        The minimum insert size for valid paired-end alignments. E.g. if -I 60 is specified and
                          a paired-end alignment consists of two 20-bp alignments in the appropriate orientation
                          with a 20-bp gap between them, that alignment is considered valid (as long as -X is also
                          satisfied). A 19-bp gap would not be valid in that case. Default: 0.

  -X/--maxins <int>        The maximum insert size for valid paired-end alignments. E.g. if -X 100 is specified and
                          a paired-end alignment consists of two 20-bp alignments in the proper orientation with a
                          60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied).
                          A 61-bp gap would not be valid in that case. Default: 500.

  --parallel <int>         (May also be --multicore <int>) Sets the number of parallel instances of Bismark to be run concurrently.
                          This forks the Bismark alignment step very early on so that each individual Spawn of Bismark processes
                          only every n-th sequence (n being set by --parallel). Once all processes have completed,
                          the individual BAM files, mapping reports, unmapped or ambiguous FastQ files are merged
                          into single files in very much the same way as they would have been generated running Bismark
                          conventionally with only a single instance.

                          If system resources are plentiful this is a viable option to speed up the alignment process
                          (we observed a near linear speed increase for up to --parallel 8 tested). However, please note
                          that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
                          Bowtie/Bowtie2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
                          and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
                          will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
                          e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
                          reduce the alignment time to ~25-30%. You have been warned.



  Bowtie 1 Reporting:

  -k <2>                   Due to the way Bismark works Bowtie will report up to 2 valid alignments. This option
                          will be used by default.

  --best                   Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum
                          (i.e. number of mismatches, or mismatches in the seed in the case if -n mode) and in
                          terms of the quality; e.g. a 1-mismatch alignment where the mismatch position has Phred
                          quality 40 is preferred over a 2-mismatch alignment where the mismatched positions both
                          have Phred quality 10. When --best is not specified, Bowtie may report alignments that
                          are sub-optimal in terms of stratum and/or quality (though an effort is made to report
                          the best alignment). --best mode also removes all strand bias. Note that --best does not
                          affect which alignments are considered "valid" by Bowtie, only which valid alignments
                          are reported by Bowtie. Bowtie is about 1-2.5 times slower when --best is specified.
                          Default: on.

  --no_best                Disables the --best option which is on by default. This can speed up the alignment process,
                          e.g. for testing purposes, but for credible results it is not recommended to disable --best.


  Output:

  --non_directional        The sequencing library was constructed in a non strand-specific manner, alignments to all four
                          bisulfite strands will be reported. Default: OFF.

                          (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary
                          to the original strands are merely theoretical and should not exist in reality. Specifying directional
                          alignments (which is the default) will only run 2 alignment threads to the original top (OT)
                          or bottom (OB) strands in parallel and report these alignments. This is the recommended option
                          for sprand-specific libraries).

  --pbat                   This options may be used for PBAT-Seq libraries (Post-Bisulfite Adapter Tagging; Kobayashi et al.,
                          PLoS Genetics, 2012). This is essentially the exact opposite of alignments in 'directional' mode,
                          as it will only launch two alignment threads to the CTOT and CTOB strands instead of the normal OT
                          and OB ones. Use this option only if you are certain that your libraries were constructed following
                          a PBAT protocol (if you don't know what PBAT-Seq is you should not specify this option). The option
                          --pbat works only for FastQ files (in both Bowtie and Bowtie 2 mode) and using uncompressed
                          temporary files only).

  --sam-no-hd              Suppress SAM header lines (starting with @). This might be useful when very large input files are
                          split up into several smaller files to run concurrently and the output files are to be merged.

  --rg_tag                 Write out a Read Group tag to the resulting SAM/BAM file. This will write the following line to the
                          SAM header: @RG PL: ILLUMINA ID:SAMPLE SM:SAMPLE ; to set ID and SM see --rg_id and --rg_sample.
                          In addition each read receives an RG:Z:RG-ID tag. Default: OFF.

  --rg_id <string>         Sets the ID field in the @RG header line. The default is 'SAMPLE'.

  --rg_sample <string>     Sets the SM field in the @RG header line; can't be set without setting --rg_id as well. The default is
                          'SAMPLE'.

  --quiet                  Print nothing besides alignments.

  --vanilla                Performs bisulfite mapping with Bowtie 1 and prints the 'old' output (as in Bismark 0.5.X) instead
                          of SAM format output.

  -un/--unmapped           Write all reads that could not be aligned to a file in the output directory. Written reads will
                          appear as they did in the input, without any translation of quality values that may have
                          taken place within Bowtie or Bismark. Paired-end reads will be written to two parallel files with _1
                          and _2 inserted in their filenames, i.e. _unmapped_reads_1.txt and unmapped_reads_2.txt. Reads
                          with more than one valid alignment with the same number of lowest mismatches (ambiguous mapping)
                          are also written to _unmapped_reads.txt unless the option --ambiguous is specified as well.

  --ambiguous              Write all reads which produce more than one valid alignment with the same number of lowest
                          mismatches or other reads that fail to align uniquely to a file in the output directory.
                          Written reads will appear as. they did in the input, without any of the translation of quality
                          values that may have taken place within Bowtie or Bismark. Paired-end reads will be written to two
                          parallel files with _1 and _2 inserted in their filenames, i.e. _ambiguous_reads_1.txt and
                          _ambiguous_reads_2.txt. These reads are not written to the file specified with --un.

  -o/--output_dir <dir>    Write all output files into this directory. By default the output files will be written into
                          the same folder as the input file(s). If the specified folder does not exist, Bismark will attempt
                          to create it first. The path to the output folder can be either relative or absolute.

  --temp_dir <dir>         Write temporary files to this directory instead of into the same directory as the input files. If
                          the specified folder does not exist, Bismark will attempt to create it first. The path to the
                          temporary folder can be either relative or absolute.

  --non_bs_mm              Optionally outputs an extra column specifying the number of non-bisulfite mismatches a read during the
                          alignment step. This option is only available for SAM format. In Bowtie 2 context, this value is
                          just the number of actual non-bisulfite mismatches and ignores potential insertions or deletions.
                          The format for single-end reads and read 1 of paired-end reads is 'XA:Z:number of mismatches'
                          and 'XB:Z:number of mismatches' for read 2 of paired-end reads.

  --gzip                   Temporary bisulfite conversion files will be written out in a GZIP compressed form to save disk
                          space. This option is available for most alignment modes but is not available for paired-end FastA
                          files. This option might be somewhat slower than writing out uncompressed files, but this awaits
                          further testing.

  --sam                    The output will be written out in SAM format instead of the default BAM format. Bismark will
                          attempt to use the path to Samtools that was specified with '--samtools_path', or, if it hasn't
                          been specified, attempt to find Samtools in the PATH. If no installation of Samtools can be found,
                          the SAM output will be compressed with GZIP instead (yielding a .sam.gz output file).

  --cram                   Writes the output to a CRAM file instead of BAM. This requires the use of Samtools 1.2 or higher.

  --cram_ref <ref_file>    CRAM output requires you to specify a reference genome as a single FastA file. If this single-FastA
                          reference file is not supplied explicitly it will be regenerated from the genome .fa sequence(s)
                          used for the Bismark run and written to a file called 'Bismark_genome_CRAM_reference.mfa' into the
                          oputput directory.

  --samtools_path          The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to be specified
                          explicitly if Samtools is in the PATH already.

  --prefix <prefix>        Prefixes <prefix> to the output filenames. Trailing dots will be replaced by a single one. For
                          example, '--prefix test' with 'file.fq' would result in the output file 'test.file.fq_bismark.sam' etc.

  -B/--basename <basename> Write all output to files starting with this base file name. For example, '--basename foo'
                          would result in the files 'foo.bam' and 'foo_SE_report.txt' (or its paired-end equivalent). Takes
                          precedence over --prefix.

  --old_flag               Only in paired-end SAM mode, uses the FLAG values used by Bismark v0.8.2 and before. In addition,
                          this options appends /1 and /2 to the read IDs for reads 1 and 2 relative to the input file. Since
                          both the appended read IDs and custom FLAG values may cause problems with some downstream tools
                          such as Picard, new defaults were implemented as of version 0.8.3.


                                              default                         old_flag
                                        ===================              ===================
                                        Read 1       Read 2              Read 1       Read 2

                                OT:         99          147                  67          131

                                OB:         83          163                 115          179

                                CTOT:      147           99                  67          131

                                CTOB:      163           83                 115          179

  --ambig_bam              For reads that have multiple alignments a random alignment is written out to a special file ending in
                          '.ambiguous.bam'. The alignments are in Bowtie2 format and do not any contain Bismark specific
                          entries such as the methylation call etc. These ambiguous BAM files are intended to be used as
                          coverage estimators for variant callers.

  --nucleotide_coverage    Calculates the mono- and di-nucleotide sequence composition of covered positions in the analysed BAM
                          file and compares it to the genomic average composition once alignments are complete by calling 'bam2nuc'.
                          Since this calculation may take a while, bam2nuc attempts to write the genomic sequence composition
                          into a file called 'genomic_nucleotide_frequencies.txt' indside the reference genome folder so it can
                          be re-used the next time round instead of calculating it once again. If a file 'nucleotide_stats.txt' is
                          found with the Bismark reports it will be automatically detected and used for the Bismark HTML report.
                          This option works only for BAM or CRAM files.


  Other:

  -h/--help                Displays this help file.

  -v/--version             Displays version information.


  BOWTIE 2 SPECIFIC OPTIONS

  --bowtie1                Uses Bowtie 1 instead of Bowtie 2, which might be a good choice for faster and very short
                          alignments. Bismark assumes that raw sequence data is adapter and/or quality trimmed where
                          appropriate. Default: off.

  --bowtie2                Default: ON. Uses Bowtie 2 instead of Bowtie 1. Bismark limits Bowtie 2 to only perform end-to-end
                          alignments, i.e. searches for alignments involving all read characters (also called
                          untrimmed or unclipped alignments). Bismark assumes that raw sequence data is adapter
                          and/or quality trimmed where appropriate. Both small (.bt2) and large (.bt2l) Bowtie 2
                          indexes are supported.

  Bowtie 2 alignment options:

  -N <int>                 Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
                          Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower)
                          but increases sensitivity. Default: 0. This option is only available for Bowtie 2 (for
                          Bowtie 1 see -n).

  -L <int>                 Sets the length of the seed substrings to align during multiseed alignment. Smaller values
                          make alignment slower but more senstive. Default: the --sensitive preset of Bowtie 2 is
                          used by default, which sets -L to 20. maximum of L can be set to 32. The length of the seed
                          would effect the alignment speed dramatically while the larger L, the faster the aligment.
                          This option is only available for Bowtie 2 (for Bowtie 1 see -l).

  --ignore-quals           When calculating a mismatch penalty, always consider the quality value at the mismatched
                          position to be the highest possible, regardless of the actual value. I.e. input is treated
                          as though all quality values are high. This is also the default behavior when the input
                          doesn't specify quality values (e.g. in -f mode). This option is invariable and on by default.


  Bowtie 2 paired-end options:

  --no-mixed               This option disables Bowtie 2's behavior to try to find alignments for the individual mates if
                          it cannot find a concordant or discordant alignment for a pair. This option is invariable and
                          and on by default.

  --no-discordant          Normally, Bowtie 2 looks for discordant alignments if it cannot find any concordant alignments.
                          A discordant alignment is an alignment where both mates align uniquely, but that does not
                          satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior
                          and it is on by default.

  --no_dovetail            It is possible, though unusual, for the mates to "dovetail", with the mates seemingly extending
                          "past" each other as in this example:

                          Mate 1:                 GTCAGCTACGATATTGTTTGGGGTGACACATTACGC
                          Mate 2:            TATGAGTCAGCTACGATATTGTTTGGGGTGACACAT
                          Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

                          By default, dovetailing is considered inconsistent with concordant alignment, but by default
                          Bismark calls Bowtie 2 with --dovetail, causing it to consider dovetailing alignments as
                          concordant. This becomes relevant whenever reads are clipped from their 5' end prior to mapping,
                          e.g. because of quality or bias issues.

                          Specify --no_dovetail to turn off this behaviour for paired-end libraries. Default: OFF.


  Bowtie 2 effort options:

  -D <int>                 Up to <int> consecutive seed extension attempts can "fail" before Bowtie 2 moves on, using
                          the alignments found so far. A seed extension "fails" if it does not yield a new best or a
                          new second-best alignment. Default: 15.

  -R <int>                 <int> is the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
                          When "re-seeding," Bowtie 2 simply chooses a new set of reads (same length, same number of
                          mismatches allowed) at different offsets and searches for more alignments. A read is considered
                          to have repetitive seeds if the total number of seed hits divided by the number of seeds
                          that aligned at least once is greater than 300. Default: 2.

  Bowtie 2 parallelization options:


  -p NTHREADS              Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores
                          and synchronize when parsing reads and outputting alignments. Searching for alignments is highly
                          parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint.
                          E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint
                          by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads
                          library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time). In addition, this option will
                          automatically use the option '--reorder', which guarantees that output SAM records are printed in
                          an order corresponding to the order of the reads in the original input file, even when -p is set
                          greater than 1 (Bismark requires the Bowtie 2 output to be this way). Specifying --reorder and
                          setting -p greater than 1 causes Bowtie 2 to run somewhat slower and use somewhat more memory then
                          if --reorder were not specified. Has no effect if -p is set to 1, since output order will naturally
                          correspond to input order in that case.

  Bowtie 2 Scoring options:

  --score_min <func>       Sets a function governing the minimum alignment score needed for an alignment to be considered
                          "valid" (i.e. good enough to report). This is a function of read length. For instance, specifying
                          L,0,-0.2 sets the minimum-score function f to f(x) = 0 + -0.2 * x, where x is the read length.
                          See also: setting function options at http://bowtie-bio.sourceforge.net/bowtie2. The default is
                          L,0,-0.2.

  --rdg <int1>,<int2>      Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets a penalty
                          of <int1> + N * <int2>. Default: 5, 3.

  --rfg <int1>,<int2>      Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of length N gets
                          a penalty of <int1> + N * <int2>. Default: 5, 3.


  Bowtie 2 Reporting options:

  -most_valid_alignments <int> This used to be the Bowtie 2 parameter -M. As of Bowtie 2 version 2.0.0 beta7 the option -M is
                          deprecated. It will be removed in subsequent versions. What used to be called -M mode is still the
                          default mode, but adjusting the -M setting is deprecated.  Use the -D and -R options to adjust the
                          effort expended to find valid alignments.

                          For reference, this used to be the old (now deprecated) description of -M:
                          Bowtie 2 searches for at most <int>+1 distinct, valid alignments for each read. The search terminates when it
                          can't find more distinct valid alignments, or when it finds <int>+1 distinct alignments, whichever
                          happens first. Only the best alignment is reported. Information from the other alignments is used to
                          estimate mapping quality and to set SAM optional fields, such as AS:i and XS:i. Increasing -M makes
                          Bowtie 2 slower, but increases the likelihood that it will pick the correct alignment for a read that
                          aligns many places. For reads that have more than <int>+1 distinct, valid alignments, Bowtie 2 does not
                          guarantee that the alignment reported is the best possible in terms of alignment score. -M is
                          always used and its default value is set to 10.


  'VANILLA' Bismark  OUTPUT:

  Single-end output format (tab-separated):

  (1) <seq-ID>
  (2) <read alignment strand>
  (3) <chromosome>
  (4) <start position>
  (5) <end position>
  (6) <observed bisulfite sequence>
  (7) <equivalent genomic sequence>
  (8) <methylation call>
  (9) <read conversion
  (10) <genome conversion>
  (11) <read quality score (Phred33)>


  Paired-end output format (tab-separated):
  (1) <seq-ID>
  (2) <read 1 alignment strand>
  (3) <chromosome>
  (4) <start position>
  (5) <end position>
  (6) <observed bisulfite sequence 1>
  (7) <equivalent genomic sequence 1>
  (8) <methylation call 1>
  (9) <observed bisulfite sequence 2>
  (10) <equivalent genomic sequence 2>
  (11) <methylation call 2>
  (12) <read 1 conversion
  (13) <genome conversion>
  (14) <read 1 quality score (Phred33)>
  (15) <read 2 quality score (Phred33)>


  Bismark SAM OUTPUT (default):

  (1) QNAME  (seq-ID)
  (2) FLAG   (this flag tries to take the strand a bisulfite read originated from into account (this is different from ordinary DNA alignment flags!))
  (3) RNAME  (chromosome)
  (4) POS    (start position)
  (5) MAPQ   (always 255 for use with Bowtie)
  (6) CIGAR
  (7) RNEXT
  (8) PNEXT
  (9) TLEN
  (10) SEQ
  (11) QUAL   (Phred33 scale)
  (12) NM-tag (edit distance to the reference)
  (13) MD-tag (base-by-base mismatches to the reference (handles indels)
  (14) XM-tag (methylation call string)
  (15) XR-tag (read conversion state for the alignment)
  (16) XG-tag (genome conversion state for the alignment)
  (17) XA/XB-tag (non-bisulfite mismatches) (optional!)

  Each read of paired-end alignments is written out in a separate line in the above format.
