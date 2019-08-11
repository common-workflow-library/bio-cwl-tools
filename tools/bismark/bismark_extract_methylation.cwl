#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2


inputs:

  genome_folder:
    type: Directory
    label: "Genome folder"
    doc: |
      "Genome folder with FASTA (fa, fasta) files.
       Bismark generated indices folder can be used also"
    inputBinding:
      position: 2
      prefix: "--genome_folder"

  bam_file:
    type: File
    label: "BAM alignment file"
    doc: "Bismark generated BAM alignment file"
    inputBinding:
      position: 3

  processes:
    type: int?
    label: "Number of Bismark instances to run"
    doc: |
      "Set the number of parallel Bismark instances to run concurrently.
       Each Bismark instance simultainously runs the methylation extractor,
       samtools stream and GZIP streams"
    inputBinding:
      position: 1
      prefix: "--multicore"


outputs:

  chg_context_file:
    type: File
    label: "CHG methylation call"
    doc: "CHG methylation call"
    outputBinding:
      glob: "CHG_context*"

  chh_context_file:
    type: File
    label: "CHH methylation call"
    doc: "CHH methylation call"
    outputBinding:
      glob: "CHH_context*"

  cpg_context_file:
    type: File
    label: "CpG methylation call"
    doc: "CpG methylation call"
    outputBinding:
      glob: "CpG_context*"

  mbias_plot:
    type: File
    label: "Methylation bias plot"
    doc: "QC data showing methylation bias across read lengths"
    outputBinding:
      glob: "*M-bias.txt"

  mbias_plot_png:
    type: File
    label: "Methylation bias plot (PNG)"
    doc: "QC data showing methylation bias across read lengths"
    outputBinding:
      glob: "*.png"

  bedgraph_coverage_file:
    type: File
    label: "Methylation statuses bedGraph coverage file"
    doc: "Coverage text file summarising cytosine methylation values in bedGraph format (tab-delimited; 0-based start coords, 1-based end coords)"
    outputBinding:
      glob: "*bedGraph.gz"

  bismark_coverage_file:
    type: File
    label: "Methylation statuses Bismark coverage file"
    doc: "Coverage text file summarising cytosine methylation values in Bismark format (tab-delimited, 1-based genomic coords)"
    outputBinding:
      glob: "*bismark.cov.gz"

  genome_wide_methylation_report:
    type: File
    label: "Genome-wide cytosine methylation report"
    doc: "Genome-wide methylation report for all cytosines in the genome"
    outputBinding:
      glob: "*CpG_report.txt"

  splitting_report:
    type: File
    label: "Methylation extraction log"
    doc: "Log file giving summary statistics about methylation extraction"
    outputBinding:
      glob: "*splitting_report.txt"


baseCommand: ["bismark_methylation_extractor", "--comprehensive", "--bedgraph", "--cytosine_report"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bismark_metadata.yaml

s:name: "bismark_extract_methylation"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bismark/bismark_extract_methylation.cwl
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
  bismark_methylation_extractor script operates on Bismark result files and extracts the methylation call
  for every single C analysed. The position of every single C will be written out to a new output file,
  depending on its context (CpG, CHG or CHH), whereby methylated Cs will be labelled as forward reads (+),
  non-methylated Cs as reverse reads (-).

  Options used:
  --comprehensive
    If strand-specific methylation is not of interest, all available methylation information can be pooled
    into a single context-dependent file (information from any of the four strands will be pooled). This
    will default to three output files (CpG-context, CHG-context and CHH-context).
  --bedgraph
    The Bismark methylation extractor can optionally also output a file in bedGraph format which uses 0-based
    genomic start and 1- based end coordinates. It will be sorted by chromosomal coordinates and looks like this:
    <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
  --genome_folder
    Bismark methylation extractor can also output a genome-wide cytosine methylation report. It is also sorted by
    chromosomal coordinates but also contains the sequence context and is in the following format:
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
    The main difference to the bedGraph or coverage output is that every cytosine on both the top and bottom strands
    will be considered irrespective of whether they were actually covered by any reads in the experiment or not.
    For this to work one has to also specify the genome that was used for the Bismark alignments using the
    option --genome_folder <path>. As for the bedGraph mode, this will only consider cytosines in CpG context.
  --cytosine_report
    After the conversion to bedGraph has completed, the option '--cytosine_report' produces a
    genome-wide methylation report for all cytosines in the genome. By default, the output uses 1-based
    chromosome coordinates (zero-based start coords are optional) and reports CpG context only (all
    cytosine context is optional). The output considers all Cs on both forward and reverse strands and
    reports their position, strand, trinucleotide content and methylation state (counts are 0 if not
    covered). The cytosine report conversion step is performed by the external module
    'coverage2cytosine'; this script needs to reside in the same folder as the bismark_methylation_extractor
    itself.

s:about: |
  DESCRIPTION

  The following is a brief description of all options to control the Bismark
  methylation extractor. The script reads in a bisulfite read alignment results file 
  produced by the Bismark bisulfite mapper (in BAM/CRAM/SAM format) and extracts the
  methylation information for individual cytosines. This information is found in the
  methylation call field which can contain the following characters:

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ~~~   X   for methylated C in CHG context                      ~~~
        ~~~   x   for not methylated C CHG                             ~~~
        ~~~   H   for methylated C in CHH context                      ~~~
        ~~~   h   for not methylated C in CHH context                  ~~~
        ~~~   Z   for methylated C in CpG context                      ~~~
        ~~~   z   for not methylated C in CpG context                  ~~~
        ~~~   U   for methylated C in Unknown context (CN or CHN       ~~~
        ~~~   u   for not methylated C in Unknown context (CN or CHN)  ~~~
        ~~~   .   for any bases not involving cytosines                ~~~
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  The methylation extractor outputs result files for cytosines in CpG, CHG and CHH
  context (this distinction is actually already made in Bismark itself). As the methylation
  information for every C analysed can produce files which easily have tens or even hundreds of
  millions of lines, file sizes can become very large and more difficult to handle. The C
  methylation info additionally splits cytosine methylation calls up into one of the four possible
  strands a given bisulfite read aligned against:

              OT      original top strand
              CTOT    complementary to original top strand

              OB      original bottom strand
              CTOB    complementary to original bottom strand

  Thus, by default twelve individual output files are being generated per input file (unless
  --comprehensive is specified, see below). The output files can be imported into a genome
  viewer, such as SeqMonk, and re-combined into a single data group if desired (in fact
  unless the bisulfite reads were generated preserving directionality it doesn't make any
  sense to look at the data in a strand-specific manner). Strand-specific output files can
  optionally be skipped, in which case only three output files for CpG, CHG or CHH context
  will be generated. For both the strand-specific and comprehensive outputs there is also
  the option to merge both non-CpG contexts (CHG and CHH) into one single non-CpG context.


  The output files are in the following format (tab delimited):

  <sequence_id>     <strand>      <chromosome>     <position>     <methylation call>


  USAGE: bismark_methylation_extractor [options] <filenames>


  ARGUMENTS:
  ==========

  <filenames>              A space-separated list of Bismark result files in SAM format from
                          which methylation information is extracted for every cytosine in
                          the reads. For alignment files in the older custom Bismark output
                          see option '--vanilla'.

  OPTIONS:

  -s/--single-end          Input file(s) are Bismark result file(s) generated from single-end
                          read data. If neither -s nor -p is set the type of experiment will
                          be determined automatically.

  -p/--paired-end          Input file(s) are Bismark result file(s) generated from paired-end
                          read data. If neither -s nor -p is set the type of experiment will
                          be determined automatically.

  --vanilla                The Bismark result input file(s) are in the old custom Bismark format
                          (up to version 0.5.x) and not in SAM format which is the default as
                          of Bismark version 0.6.x or higher. Default: OFF.

  --no_overlap             For paired-end reads it is theoretically possible that read_1 and
                          read_2 overlap. This option avoids scoring overlapping methylation
                          calls twice (only methylation calls of read 1 are used for in the process
                          since read 1 has historically higher quality basecalls than read 2).
                          Whilst this option removes a bias towards more methylation calls
                          in the center of sequenced fragments it may de facto remove a sizable
                          proportion of the data. This option is on by default for paired-end data
                          but can be disabled using --include_overlap. Default: ON.

  --include_overlap        For paired-end data all methylation calls will be extracted irrespective of
                          of whether they overlap or not. Default: OFF.

  --ignore <int>           Ignore the first <int> bp from the 5' end of Read 1 (or single-end alignment
                          files) when processing the methylation call string. This can remove e.g. a
                          restriction enzyme site at the start of each read or any other source of
                          bias (such as PBAT-Seq data).

  --ignore_r2 <int>        Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing
                          results only. Since the first couple of bases in Read 2 of BS-Seq experiments
                          show a severe bias towards non-methylation as a result of end-repairing
                          sonicated fragments with unmethylated cytosines (see M-bias plot), it is
                          recommended that the first couple of bp of Read 2 are removed before
                          starting downstream analysis. Please see the section on M-bias plots in the
                          Bismark User Guide for more details.

  --ignore_3prime <int>    Ignore the last <int> bp from the 3' end of Read 1 (or single-end alignment
                          files) when processing the methylation call string. This can remove unwanted
                          biases from the end of reads.

  --ignore_3prime_r2 <int> Ignore the last <int> bp from the 3' end of Read 2 of paired-end sequencing
                          results only. This can remove unwanted biases from the end of reads.

  --comprehensive          Specifying this option will merge all four possible strand-specific
                          methylation info into context-dependent output files. The default

                          contexts are:
                            - CpG context
                            - CHG context
                            - CHH context

  --merge_non_CpG          This will produce two output files (in --comprehensive mode) or eight
                          strand-specific output files (default) for Cs in
                            - CpG context
                            - non-CpG context

  --report                 Prints out a short methylation summary as well as the paramaters used to run
                          this script. Default: ON.

  --no_header              Suppresses the Bismark version header line in all output files for more convenient
                          batch processing.

  -o/--output DIR          Allows specification of a different output directory (absolute or relative
                          path). If not specified explicitly, the output will be written to the current directory.

  --samtools_path          The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to be specified
                          explicitly if Samtools is in the PATH already.

  --gzip                   The methylation extractor files (CpG_OT_..., CpG_OB_... etc) will be written out in
                          a GZIP compressed form to save disk space. This option is also passed on to the genome-wide
                          cytosine report. BedGraph and coverage files will be written out as .gz by default.

  --version                Displays version information.

  -h/--help                Displays this help file and exits.

  --mbias_only             The methylation extractor will read the entire file but only output the M-bias table and plots as 
                          well as a report (optional) and then quit. Default: OFF.

  --mbias_off              The methylation extractor will process the entire file as usual but doesn't write out any M-bias report.
                          Only recommended for users who deliberately want to keep an earlier version of the M-bias report. 
                          Default: OFF.

  --parallel <int>         May also be --multicore <int>. Sets the number of cores to be used for the methylation extraction process.
                          If system resources are plentiful this is a viable option to speed up the extraction process (we observed a
                          near linear speed increase for up to 10 cores used). Please note that a typical process of extracting a BAM file
                          and writing out '.gz' output streams will in fact use ~3 cores per value of --parallel <int>
                          specified (1 for the methylation extractor itself, 1 for a Samtools stream, 1 for GZIP stream), so
                          --parallel 10 is likely to use around 30 cores of system resources. This option has no bearing
                          on the bismark2bedGraph or genome-wide cytosine report processes.

  --yacht                  This option (Yet Another Context Hunting Tool) writes out additional information about the read a methylation
                          call belongs to, and its output is meant to be fed into the NOMe_filtering script. This option writes out
                          a single 'any_C_context' file that contains all methylation calls for a read consecutively. Its intended use
                          is single-cell NOMe-Seq data, and thus this option works only in single-end mode (paired-end reads often suffer
                          from chimaera problems...)

                          --yacht will add three additional columns to the standard methylation call files:
                                      <read start>  <read end>  <read orientation>
                          For forward reads (+ orientation) the start position is the left-most position wheras for reverse reads
                          (- orientation) it is the rightmost position.


  bedGraph specific options:
  ==========================

  --bedGraph               After finishing the methylation extraction, the methylation output is written into a
                          sorted bedGraph file that reports the position of a given cytosine and its methylation 
                          state (in %, see details below). The methylation extractor output is temporarily split up into
                          temporary files, one per chromosome (written into the current directory or folder
                          specified with -o/--output); these temp files are then used for sorting and deleted
                          afterwards. By default, only cytosines in CpG context will be sorted. The option
                          '--CX_context' may be used to report all cytosines irrespective of sequence context
                          (this will take MUCH longer!). The default folder for temporary files during the sorting
                          process is the output directory. The bedGraph conversion step is performed by the external
                          module 'bismark2bedGraph'; this script needs to reside in the same folder as the 
                          bismark_methylation_extractor itself.

  --zero_based             Write out an additional coverage file (ending in .zero.cov) that uses 0-based genomic start
                          and 1-based genomic end coordinates (zero-based, half-open), like used in the bedGraph file,
                          instead of using 1-based coordinates throughout. Default: OFF.


  --cutoff [threshold]     The minimum number of times any methylation state (methylated or unmethylated) has to be seen
                          for a nucleotide before its methylation percentage is reported. Default: 1.

  --remove_spaces          Replaces whitespaces in the sequence ID field with underscores to allow sorting.

  --CX/--CX_context        The sorted bedGraph output file contains information on every single cytosine that was covered
                          in the experiment irrespective of its sequence context. This applies to both forward and
                          reverse strands. Please be aware that this option may generate large temporary and output files
                          and may take a long time to sort (up to many hours). Default: OFF.
                          (i.e. Default = CpG context only).

  --buffer_size <string>   This allows you to specify the main memory sort buffer when sorting the methylation information.
                          Either specify a percentage of physical memory by appending % (e.g. --buffer_size 50%) or
        a multiple of 1024 bytes, e.g. 'K' multiplies by 1024, 'M' by 1048576 and so on for 'T' etc. 
                          (e.g. --buffer_size 20G). For more information on sort type 'info sort' on a command line.
                          Defaults to 2G.

  --scaffolds/--gazillion  Users working with unfinished genomes sporting tens or even hundreds of thousands of
                          scaffolds/contigs/chromosomes frequently encountered errors with pre-sorting reads to
                          individual chromosome files. These errors were caused by the operating system's limit
                          of the number of filehandle that can be written to at any one time (typically 1024; to
                          find out this limit on Linux, type: ulimit -a).
                          To bypass the limitation of open filehandles, the option --scaffolds does not pre-sort
                          methylation calls into individual chromosome files. Instead, all input files are
                          temporarily merged into a single file (unless there is only a single file), and this
                          file will then be sorted by both chromosome AND position using the Unix sort command.
                          Please be aware that this option might take a looooong time to complete, depending on
                          the size of the input files, and the memory you allocate to this process (see --buffer_size).
                          Nevertheless, it seems to be working.

  --ample_memory           Using this option will not sort chromosomal positions using the UNIX 'sort' command, but will
                          instead use two arrays to sort methylated and unmethylated calls. This may result in a faster
                          sorting process of very large files, but this comes at the cost of a larger memory footprint
                          (two arrays of the length of the largest human chromosome 1 (~250M bp) consume around 16GB
                          of RAM). Due to overheads in creating and looping through these arrays it seems that it will
                          actually be *slower* for small files (few million alignments), and we are currently testing at
                          which point it is advisable to use this option. Note that --ample_memory is not compatible
                          with options '--scaffolds/--gazillion' (as it requires pre-sorted files to begin with).



  Genome-wide cytosine methylation report specific options:
  =========================================================

  --cytosine_report        After the conversion to bedGraph has completed, the option '--cytosine_report' produces a
                          genome-wide methylation report for all cytosines in the genome. By default, the output uses 1-based
                          chromosome coordinates (zero-based start coords are optional) and reports CpG context only (all
                          cytosine context is optional). The output considers all Cs on both forward and reverse strands and
                          reports their position, strand, trinucleotide content and methylation state (counts are 0 if not
                          covered). The cytosine report conversion step is performed by the external module
                          'coverage2cytosine'; this script needs to reside in the same folder as the bismark_methylation_extractor
                          itself.

  --CX/--CX_context        The output file contains information on every single cytosine in the genome irrespective of
                          its context. This applies to both forward and reverse strands. Please be aware that this will
                          generate output files with > 1.1 billion lines for a mammalian genome such as human or mouse.
                          Default: OFF (i.e. Default = CpG context only).

  --zero_based             Uses 0-based genomic coordinates instead of 1-based coordinates. Default: OFF.

  --genome_folder <path>   Enter the genome folder you wish to use to extract sequences from (full path only). Accepted
                          formats are FastA files ending with '.fa' or '.fasta'. Specifying a genome folder path is mandatory.

  --split_by_chromosome    Writes the output into individual files for each chromosome instead of a single output file. Files
                          will be named to include the input filename and the chromosome number.

  OUTPUT:

  The bismark_methylation_extractor output is in the form:
  ========================================================
  <seq-ID>  <methylation state*>  <chromosome>  <start position (= end position)>  <methylation call>

  * Methylated cytosines receive a '+' orientation,
  * Unmethylated cytosines receive a '-' orientation.

  The bismark_methylation_extractor output with --yacht (optional) specified is in the form:
  ==========================================================================================
  <seq-ID>  <methylation state*>  <chromosome>  <start position (= end position)>  <methylation call>  <read start>  <read end>  <read orientation>

  * Methylated cytosines receive a '+' orientation,
  * Unmethylated cytosines receive a '-' orientation.


  The bedGraph output (optional) looks like this (tab-delimited; 0-based start coords, 1-based end coords):
  =========================================================================================================
  track type=bedGraph (header line)
  <chromosome>  <start position>  <end position>  <methylation percentage>


  The coverage output looks like this (tab-delimited, 1-based genomic coords; zero-based half-open coordinates available with '--zero_based'):
  ============================================================================================================================================
  <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>


  The genome-wide cytosine methylation output file is tab-delimited in the following format:
  ==========================================================================================
  <chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>
