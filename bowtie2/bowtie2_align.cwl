#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        if (Array.isArray(inputs.filelist) && inputs.filelist.length > 0){
          return inputs.filelist[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sam";
        } else
          if (inputs.filelist != null){
            return inputs.filelist.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sam";
          } else
            if (Array.isArray(inputs.filelist_mates) && inputs.filelist_mates.length > 0){
              return inputs.filelist_mates[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sam";
            } else
              if (inputs.filelist_mates != null){
                return inputs.filelist_mates.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".sam";
              } else {
                return null;
              }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie2:v2.3.0
  dockerFile: >
    $import: ./dockerfiles/bowtie2-Dockerfile
- class: SoftwareRequirement
  packages:
      bowtie2:
        specs: [ "http://identifiers.org/biotools/bowtie2" ]
        version: [ "2.3.0" ]

inputs:

  indices_file:
    type: File?
    doc: File with secondaryFiles containing all the indices
    inputBinding:
      position: 81
      prefix: -x
      valueFrom: $(self.path.split('.').slice(0,-2).join('.'))

  indices_folder:
    type: Directory?
    doc: "Folder with indices files"
    inputBinding:
      position: 81
      prefix: '-x'
      valueFrom: |
        ${
            for (var i = 0; i < self.listing.length; i++) {
                if (self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.bt2' ||
                    self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.bt2l'){
                  return self.listing[i].path.split('.').slice(0,-3).join('.');
                }
            }
            return null;
        }

  filelist:
    type:
    - "null"
    - File
    - type: array
      items: File
    doc: |
      {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
      <m1>       Files with #1 mates, paired with files in <m2>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <m2>       Files with #2 mates, paired with files in <m1>.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
      <r>        Files with unpaired reads.
                 Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
    inputBinding:
      itemSeparator: ","
      position: 83

  filelist_mates:
    type:
    - "null"
    - File
    - type: array
      items: File
    inputBinding:
      itemSeparator: ","
      position: 85

  output_filename:
    type: string
    inputBinding:
      position: 90
      prefix: "-S"
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ""
    doc: |
      File for SAM output (default: stdout)

  q:
    type:
    - "null"
    - boolean
    doc: "query input files are FASTQ .fq/.fastq (default)"
    inputBinding:
      position: 1
      prefix: '-q'

  qseq:
    type:
    - "null"
    - boolean
    doc: "query input files are in Illumina's qseq format"
    inputBinding:
      position: 1
      prefix: '--qseq'

  f:
    type:
    - "null"
    - boolean
    doc: "query input files are (multi-)FASTA .fa/.mfa"
    inputBinding:
      position: 1
      prefix: '-f'

  raw:
    type:
    - "null"
    - boolean
    doc: "query input files are raw one-sequence-per-line"
    inputBinding:
      position: 1
      prefix: '-r'

  c:
    type:
    - "null"
    - boolean
    doc: "<m1>, <m2>, <r> are sequences themselves, not files"
    inputBinding:
      position: 1
      prefix: '-c'

  s:
    type:
    - "null"
    - int
    doc: |
      skip the first <int> reads/pairs in the input (none)
    inputBinding:
      position: 2
      prefix: '-s'

  u:
    type:
    - "null"
    - int
    doc: |
      stop after first <int> reads/pairs (no limit)
    inputBinding:
      position: 3
      prefix: '-u'

  clip_5p_end:
    type: int?
    doc: |
      trim <int> bases from 5'/left end of reads (0)
    inputBinding:
      position: 4
      prefix: '-5'

  clip_3p_end:
    type: int?
    doc: |
      trim <int> bases from 3'/right end of reads (0)
    inputBinding:
      position: 5
      prefix: '-3'

  phred33_quals:
    type:
    - "null"
    - boolean
    doc: "qualities are Phred+33 (default)"
    inputBinding:
      position: 6
      prefix: '--phred33'

  phred64_quals:
    type:
    - "null"
    - boolean
    doc: "qualities are Phred+64"
    inputBinding:
      position: 6
      prefix: '--phred64'

  integer_quals:
    type:
    - "null"
    - boolean
    doc: "qualities encoded as space-delimited integers"
    inputBinding:
      position: 6
      prefix: '--int-quals'

  n:
    type:
    - "null"
    - int
    doc: |
      max # mismatches in seed alignment; can be 0 or 1 (0)
    inputBinding:
      position: 7
      prefix: '-N'

  l:
    type:
    - "null"
    - int
    doc: |
      length of seed substrings; must be >3, <32 (22)
    inputBinding:
      position: 8
      prefix: '-L'

  i:
    type:
    - "null"
    - int
    doc: |
      interval between seed substrings w/r/t read len (S,1,1.15)
    inputBinding:
      position: 9
      prefix: '-i'

  n_ceil:
    type:
    - "null"
    - string
    doc: >
      func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
    inputBinding:
      position: 10
      prefix: '--n-ceil'

  dpad:
    type:
    - "null"
    - int
    doc: |
      include <int> extra ref chars on sides of DP table (15)
    inputBinding:
      position: 11
      prefix: '--dpad'

  gbar:
    type:
    - "null"
    - int
    doc: |
      disallow gaps within <int> nucs of read extremes (4)
    inputBinding:
      position: 12
      prefix: '--gbar'

  ignore_quals:
    type:
    - "null"
    - boolean
    doc: |
      treat all quality values as 30 on Phred scale (off)
    inputBinding:
      position: 13
      prefix: '--ignore-quals'

  nofw:
    type:
    - "null"
    - boolean
    doc: |
      do not align forward (original) version of read (off)
    inputBinding:
      position: 14
      prefix: '--nofw'

  norc:
    type:
    - "null"
    - boolean
    doc: |
      do not align reverse-complement version of read (off)
    inputBinding:
      position: 15
      prefix: '--norc'

  no_1mm_upfront:
    type:
    - "null"
    - boolean
    doc: |
      do not allow 1 mismatch alignments before attempting to scan for the optimal seeded alignments
    inputBinding:
      position: 16
      prefix: '--no-1mm-upfront'

  end_to_end:
    type:
    - "null"
    - boolean
    doc: |
      entire read must align; no clipping (on)
      Options:
        --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
        --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
        --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
        --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    inputBinding:
      position: 17
      prefix: '--end-to-end'

  end_to_end_very_fast:
    type:
    - "null"
    - boolean
    doc: |
      Option for end_to_end:
        --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
    inputBinding:
      position: 18
      prefix: '--very-fast'

  end_to_end_fast:
    type:
    - "null"
    - boolean
    doc: |
      Option for end_to_end:
        --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
    inputBinding:
      position: 18
      prefix: '--fast'

  end_to_end_sensitive:
    type:
    - "null"
    - boolean
    doc: |
      Option for end_to_end:
        --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)
    inputBinding:
      position: 18
      prefix: '--sensitive'

  end_to_end_very_sensitive:
    type:
    - "null"
    - boolean
    doc: |
      Option for end_to_end:
        --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    inputBinding:
      position: 18
      prefix: '--very-sensitive'

  local:
    type:
    - "null"
    - boolean
    doc: |
      local alignment; ends might be soft clipped (off)
      Options:
        --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
        --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
        --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
        --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    inputBinding:
      position: 19
      prefix: '--local'

  local_very_fast_local:
    type:
    - "null"
    - boolean
    doc: |
      Option for local:
        --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
    inputBinding:
      position: 20
      prefix: '--very-fast-local'

  local_fast_local:
    type:
    - "null"
    - boolean
    doc: |
      Option for local:
        --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75
    inputBinding:
      position: 20
      prefix: '--fast-local'

  local_sensitive_local:
    type:
    - "null"
    - boolean
    doc: |
      Option for local:
        --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
    inputBinding:
      position: 20
      prefix: '--sensitive-local'

  local_very_sensitive_local:
    type:
    - "null"
    - boolean
    doc: |
      Option for local:
        --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    inputBinding:
      position: 20
      prefix: '--very-sensitive-local'

  ma:
    type:
    - "null"
    - int
    doc: |
      match bonus (0 for --end-to-end, 2 for --local)
    inputBinding:
      position: 21
      prefix: '--ma'

  mp:
    type:
    - "null"
    - int
    doc: |
      max penalty for mismatch; lower qual = lower penalty (6)
    inputBinding:
      position: 22
      prefix: '--mp'

  np:
    type:
    - "null"
    - int
    doc: |
      penalty for non-A/C/G/Ts in read/ref (1)
    inputBinding:
      position: 23
      prefix: '--np'

  rdg:
    type:
    - "null"
    - int[]
    doc: |
      read gap open, extend penalties (5,3)
    inputBinding:
      position: 24
      itemSeparator: ","
      prefix: '--rdg'

  rfg:
    type:
    - "null"
    - int[]
    doc: |
      reference gap open, extend penalties (5,3)
    inputBinding:
      position: 25
      itemSeparator: ","
      prefix: '--rfg'

  score_min:
    type:
    - "null"
    - string
    doc: |
      min acceptable alignment score w/r/t read length (G,20,8 for local, L,-0.6,-0.6 for end-to-end)
    inputBinding:
      position: 26
      prefix: '--score-min'

  k:
    type:
    - "null"
    - int
    doc: |
      report up to <int> alns per read; MAPQ not meaningful
    inputBinding:
      position: 27
      prefix: '-k'

  a:
    type:
    - "null"
    - boolean
    doc: |
      report all alignments; very slow, MAPQ not meaningful
    inputBinding:
      position: 27
      prefix: '-a'

  d:
    type:
    - "null"
    - int
    doc: |
      give up extending after <int> failed extends in a row (15)
    inputBinding:
      position: 28
      prefix: '-D'

  r:
    type:
    - "null"
    - int
    doc: |
      for reads w/ repetitive seeds, try <int> sets of seeds (2)
    inputBinding:
      position: 29
      prefix: '-R'

  minins:
    type:
    - "null"
    - int
    doc: |
      minimum fragment length (0)
    inputBinding:
      position: 30
      prefix: '--minins'

  maxins:
    type:
    - "null"
    - int
    doc: |
      maxins fragment length (0)
    inputBinding:
      position: 30
      prefix: '--maxins'

  fr:
    type:
    - "null"
    - boolean
    doc: |
      -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
    inputBinding:
      position: 31
      prefix: '--fr'

  rf:
    type:
    - "null"
    - boolean
    doc: |
      -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
    inputBinding:
      position: 31
      prefix: '--rf'

  ff:
    type:
    - "null"
    - boolean
    doc: |
      -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
    inputBinding:
      position: 31
      prefix: '--ff'

  no_mixed:
    type:
    - "null"
    - boolean
    doc: |
      suppress unpaired alignments for paired reads
    inputBinding:
      position: 32
      prefix: '--no-mixed'

  no_discordant:
    type:
    - "null"
    - boolean
    doc: |
      suppress discordant alignments for paired reads
    inputBinding:
      position: 33
      prefix: '--no-discordant'

  no_dovetail:
    type:
    - "null"
    - boolean
    doc: |
      not concordant when mates extend past each other
    inputBinding:
      position: 34
      prefix: '--no-dovetail'

  no_contain:
    type:
    - "null"
    - boolean
    doc: |
      not concordant when one mate alignment contains other
    inputBinding:
      position: 35
      prefix: '--no-contain'

  no_overlap:
    type:
    - "null"
    - boolean
    doc: |
      not concordant when mates overlap at all
    inputBinding:
      position: 36
      prefix: '--no-overlap'

  t:
    type:
    - "null"
    - boolean
    doc: |
      print wall-clock time taken by search phases
    inputBinding:
      position: 37
      prefix: '-t'

  un:
    type:
    - "null"
    - string
    doc: |
      write unpaired reads that didn't align to <path>
    inputBinding:
      position: 38
      prefix: '--un'

  al:
    type:
    - "null"
    - string
    doc: |
      write unpaired reads that aligned at least once to <path>
    inputBinding:
      position: 39
      prefix: '--al'

  un_conc:
    type:
    - "null"
    - string
    doc: |
      write pairs that didn't align concordantly to <path>
    inputBinding:
      position: 40
      prefix: '--un-conc'

  al_conc:
    type:
    - "null"
    - string
    doc: |
      write pairs that aligned concordantly at least once to <path>
    inputBinding:
      position: 41
      prefix: '--al-conc'

  quiet:
    type:
    - "null"
    - boolean
    doc: "print nothing to stderr except serious errors"
    inputBinding:
      position: 42
      prefix: '--quiet'

  met_file:
    type:
    - "null"
    - string
    doc: |
      send metrics to file at <path> (off)
    inputBinding:
      position: 43
      prefix: '--met-file'

  met_stderr:
    type:
    - "null"
    - boolean
    doc: "send metrics to stderr (off)"
    inputBinding:
      position: 44
      prefix: '--met-stderr'

  met:
    type:
    - "null"
    - int
    doc: |
      report internal counters & metrics every <int> secs (1)
    inputBinding:
      position: 45
      prefix: '--met'

  no_unal:
    type:
    - "null"
    - boolean
    doc: "suppress SAM records for unaligned reads"
    inputBinding:
      position: 46
      prefix: '--no-unal'

  no_head:
    type:
    - "null"
    - boolean
    doc: "suppress header lines, i.e. lines starting with @"
    inputBinding:
      position: 47
      prefix: '--no-head'

  no_sq:
    type:
    - "null"
    - boolean
    doc: "suppress @SQ header lines"
    inputBinding:
      position: 48
      prefix: '--no-sq'

  rg_id:
    type:
    - "null"
    - string
    doc: "set read group id, reflected in @RG line and RG:Z: opt field"
    inputBinding:
      position: 49
      prefix: '--rg-id'

  rg:
    type:
    - "null"
    - string
    doc: |
      add <text> ("lab:value") to @RG line of SAM header.
      Note: @RG line only printed when --rg-id is set.
    inputBinding:
      position: 50
      prefix: '--rg'

  omit_sec_seq:
    type:
    - "null"
    - boolean
    doc: "put '*' in SEQ and QUAL fields for secondary alignments"
    inputBinding:
      position: 51
      prefix: '--omit-sec-seq'

  threads:
    type:
    - "null"
    - int
    doc: |
      number of alignment threads to launch (1)
    inputBinding:
      position: 52
      prefix: '-p'

  reorder:
    type:
    - "null"
    - boolean
    doc: |
      force SAM output order to match order of input reads
    inputBinding:
      position: 53
      prefix: '--reorder'

  mm:
    type:
    - "null"
    - boolean
    doc: "use memory-mapped I/O for index; many 'bowtie's can share"
    inputBinding:
      position: 54
      prefix: '--mm'

  qc_filter:
    type:
    - "null"
    - boolean
    doc: "filter out reads that are bad according to QSEQ filter"
    inputBinding:
      position: 55
      prefix: '--qc-filter'
#
  seed:
    type:
    - "null"
    - int
    doc: |
      seed for random number generator (0)
    inputBinding:
      position: 56
      prefix: '--seed'

  non_deterministic:
    type:
    - "null"
    - boolean
    doc: "seed rand. gen. arbitrarily instead of using read attributes"
    inputBinding:
      position: 57
      prefix: '--non-deterministic'

outputs:

  output:
    type: File
    outputBinding:
      glob: |
        ${
           if (inputs.output_filename == ""){
             return default_output_filename();
           } else {
             return inputs.output_filename;
           }
        }

  output_log:
    type: File
    outputBinding:
      glob: |
        ${
           if (inputs.output_filename == ""){
             return default_output_filename().split('.').slice(0,-1).join('.') + ".log";
           } else {
             return inputs.output_filename.split('.').slice(0,-1).join('.') + ".log";
           }
        }

baseCommand:
  - bowtie2

arguments:
  - valueFrom: |
      ${
        if (inputs.filelist && inputs.filelist_mates){
          return "-1";
        } else if (inputs.filelist){
          return "-U";
        } else {
          return null;
        }
      }
    position: 82
  - valueFrom: |
      ${
        if (inputs.filelist && inputs.filelist_mates){
          return "-2";
        } else if (inputs.filelist_mates){
          return "-U";
        } else {
          return null;
        }
      }
    position: 84
  - valueFrom: |
      ${
        if (inputs.output_filename == ""){
          return ' 2> ' + default_output_filename().split('.').slice(0,-1).join('.') + '.log';
        } else {
          return ' 2> ' + inputs.output_filename.split('.').slice(0,-1).join('.') + '.log';
        }
      }
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bowtie2_align"
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Tool runs bowtie aligner to align input FASTQ file(s) to reference genome
