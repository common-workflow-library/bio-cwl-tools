#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) {
      if (inputs.output_filename != "" && !ext){
         return inputs.output_filename;
      }
      ext = ext || ".sam";
      var root = "";
      if (inputs.output_filename != ""){
         root = inputs.output_filename.split('.').slice(0,-1).join('.');
         return (root == "")?inputs.output_filename+ext:root+ext;
      } else
        if (Array.isArray(inputs.upstream_filelist) && inputs.upstream_filelist.length > 0){
          root = inputs.upstream_filelist[0].basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.upstream_filelist[0].basename+ext:root+ext;
        } else
          if (inputs.upstream_filelist != null){
            root = inputs.upstream_filelist.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.upstream_filelist.basename+ext:root+ext;
          } else
            if (Array.isArray(inputs.downstream_filelist) && inputs.downstream_filelist.length > 0){
              root = inputs.downstream_filelist[0].basename.split('.').slice(0,-1).join('.');
              return (root == "")?inputs.downstream_filelist[0].basename+ext:root+ext;
            } else
              if (inputs.downstream_filelist != null){
                root = inputs.downstream_filelist.basename.split('.').slice(0,-1).join('.');
                return (root == "")?inputs.downstream_filelist.basename+ext:root+ext;
              } else
                if (Array.isArray(inputs.crossbow_filelist) && inputs.crossbow_filelist.length > 0){
                  root = inputs.crossbow_filelist[0].basename.split('.').slice(0,-1).join('.');
                  return (root == "")?inputs.crossbow_filelist[0].basename+ext:root+ext;
                } else
                  if (inputs.crossbow_filelist != null){
                    root = inputs.crossbow_filelist.basename.split('.').slice(0,-1).join('.');
                    return (root == "")?inputs.crossbow_filelist.basename+ext:root+ext;
                  } else {
                    return null;
                  }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bowtie:v1.2.0
- class: SoftwareRequirement
  packages:
    bowtie:
      specs: [ "http://identifiers.org/biotools/bowtie" ]
      version: [ "1.2.0" ]

inputs:

  indices_folder:
    type: Directory
    inputBinding:
      position: 81
      valueFrom: |
        ${
            for (var i = 0; i < self.listing.length; i++) {
                if (self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.ebwt' ||
                    self.listing[i].path.split('.').slice(-3).join('.') == 'rev.1.ebwtl'){
                  return self.listing[i].path.split('.').slice(0,-3).join('.');
                }
            }
            return null;
        }
    doc: |
      Folder with Bowtie indices

  upstream_filelist:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 83
    doc: |
      Comma-separated list of files containing upstream mates (or the
      sequences themselves, if -c is set) paired with mates in <m2>

  downstream_filelist:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 85
    doc: |
      Comma-separated list of files containing downstream mates (or the
      sequences themselves if -c is set) paired with mates in <m1>

  crossbow_filelist:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 86
      prefix: "-12"
    doc: |
      Comma-separated list of files containing Crossbow-style reads.
      Can be a mixture of paired and unpaired.  Specify "-"for stdin.

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 90
      valueFrom: $(default_output_filename())
    default: ""
    doc: |
      Generates default output filename on the base of upstream_filelist/downstream_filelist files

  q:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-q'
    doc: |
      query input files are FASTQ .fq/.fastq (default)

  f:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-f'
    doc: |
      query input files are (multi-)FASTA .fa/.mfa

  r:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-r'
    doc: |
      query input files are raw one-sequence-per-line

  c:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-c'
    doc: |
      query sequences given on cmd line (as <mates>, <singles>)

  C:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-C'
    doc: |
      reads and index are in colorspace

  Q:
    type:
      - "null"
      - File
    inputBinding:
      position: 1
      prefix: '-Q'
    doc: |
      QV file(s) corresponding to CSFASTA inputs; use with -f -C

  Q1:
    type:
      - "null"
      - File
    inputBinding:
      position: 1
      prefix: '--Q1'
    doc: |
      --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively

  Q2:
    type:
      - "null"
      - File
    inputBinding:
      position: 1
      prefix: '--Q2'
    doc: |
      --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively

  s:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-s'
    doc: |
      --skip <int>    skip the first <int> reads/pairs in the input

  u:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-u'
    doc: |
      --qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)

  clip_5p_end:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-5'
    doc: |
      --trim5 <int>   trim <int> bases from 5' (left) end of reads

  clip_3p_end:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-3'
    doc: |
      --trim3 <int>   trim <int> bases from 3' (right) end of reads

  phred33_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--phred33-quals'
    doc: |
      input quals are Phred+33 (default)

  phred64_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--phred64-quals'
    doc: |
      input quals are Phred+64 (same as --solexa1.3-quals)

  solexa_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--solexa-quals'
    doc: |
      input quals are from GA Pipeline ver. < 1.3

  solexa1.3_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--solexa1.3-quals'
    doc: |
      input quals are from GA Pipeline ver. >= 1.3

  integer_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--integer-quals'
    doc: |
      qualities are given as space-separated integers (not ASCII)

  large_index:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--large-index'
    doc: |
      force usage of a 'large' index, even if a small one is present

  v:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-v'
    doc: |
      report end-to-end hits w/ <=v mismatches; ignore qualities

  n:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-n'
    doc: |
      --seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)

  e:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-e'
    doc: |
      --maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)

  l:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-l'
    doc: |
      --seedlen <int> seed length for -n (default: 28)

  nomaqround:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--nomaqround'
    doc: |
      disable Maq-like quality rounding for -n (nearest 10 <= 30)

  I:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-I'
    doc: |
      --minins <int>  minimum insert size for paired-end alignment (default: 0)

  X:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-X'
    doc: |
      --maxins <int>  maximum insert size for paired-end alignment (default: 250)

  fr:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--fr'
    doc: |
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)

  rf:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--rf'
    doc: |
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)

  ff:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--ff'
    doc: |
      --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)

  nofw:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--nofw'
    doc: |
      --norc      do not align to forward/reverse-complement reference strand

  maxbts:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--maxbts'
    doc: |
      <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)

  pairtries:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--pairtries'
    doc: |
      <int>  max # attempts to find mate for anchor hit (default: 100)

  y:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-y'
    doc: |
      --tryhard try hard to find valid alignments, at the expense of speed

  chunkmbs:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--chunkmbs'
    doc: |
      <int>   max megabytes of RAM for best-first search frames (def: 64)

  reads_per_batch:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--reads-per-batch'
    doc: |
      # of reads to read from input file at once (default: 16)

  k:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-k'
    doc: |
      <int>           report up to <int> good alignments per read (default: 1)

  a:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-a'
    doc: |
      --all           report all alignments per read (much slower than low -k)

  m:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-m'
    doc: |
      <int>           suppress all alignments if > <int> exist (def: no limit)

  M:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-M'
    doc: |
      <int> like -m, but reports 1 random hit (MAPQ=0)

  best:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--best'
    doc: |
      hits guaranteed best stratum; ties broken by quality

  strata:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--strata'
    doc: |
      hits in sub-optimal strata aren't reported (requires --best)

  t:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-t'
    doc: |
      --time          print wall-clock time taken by search phases

  B:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-B'
    doc: |
      --offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)

  quiet:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--quiet'
    doc: |
      print nothing but the alignments

  refout:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--refout'
    doc: |
      write alignments to files refXXXXX.map, 1 map per reference

  refidx:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--refidx'
    doc: |
      refer to ref. seqs by 0-based index rather than name

  al:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      prefix: '--al'
    doc: |
      <fname>       write aligned reads/pairs to file(s) <fname>

  un:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      prefix: '--un'
    doc: |
      <fname>       write unaligned reads/pairs to file(s) <fname>

  max:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      prefix: '--max'
    doc: |
      <fname>      write reads/pairs over -m limit to file(s) <fname>

  suppress:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      prefix: '--suppress'
    doc: |
      <cols>  suppresses given columns (comma-delim'ed) in default output

  fullref:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--fullref'
    doc: |
      write entire ref name (default: only up to 1st space)

  snpphred:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--snpphred'
    doc: |
      <int>   Phred penalty for SNP when decoding colorspace (def: 30)

  snpfrac:
    type:
      - "null"
      - float
    inputBinding:
      position: 1
      prefix: '--snpfrac'
    doc: |
      <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred

  col_cseq:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--col-cseq'
    doc: |
      print aligned colorspace seqs as colors, not decoded bases

  col_cqual:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--col-cqual'
    doc: |
      print original colorspace quals, not decoded quals

  col_keepends:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--col-keepends'
    doc: |
      keep nucleotides at extreme ends of decoded alignment

  sam:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '-S'
    doc: |
      --sam write hits in SAM format

  mapq:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--mapq'
    doc: |
      <int>       default mapping quality (MAPQ) to print for SAM alignments

  sam_nohead:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--sam-nohead'
    doc: |
      supppress header lines (starting with @) for SAM output

  sam_nosq:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--sam-nosq'
    doc: |
      supppress @SQ header lines for SAM output

  sam_RG:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      prefix: '--sam-RG'
    doc: |
      <text>    add <text> (usually "lab=value") to @RG line of SAM header

  o:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-o'
    doc: |
      --offrate <int> override offrate of index; must be >= index's offrate

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '-p'
    doc: |
      --threads <int> number of alignment threads to launch (default: 1)

  mm:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--mm'
    doc: |
      use memory-mapped I/O for index; many 'bowtie's can share

  shmem:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--shmem'
    doc: |
      use shared mem for index

  seed:
    type:
      - "null"
      - int
    inputBinding:
      position: 1
      prefix: '--seed'
    doc: |
      <int>       seed for random number generator

  verbose:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 1
      prefix: '--verbose'
    doc: |
      verbose output (for debugging)

outputs:

  sam_file:
    type: File?
    outputBinding:
      glob: $(default_output_filename())

  refout_file:
    type: File[]?
    outputBinding:
      glob: "*ref*.map*"

  unaligned_file:
    type: File?
    outputBinding:
      glob: $(inputs.un)

  aligned_file:
    type: File?
    outputBinding:
      glob: $(inputs.al)

  multimapped_file:
    type: File?
    outputBinding:
      glob: $(inputs.max)

  log_file:
    type: File
    outputBinding:
      glob: $(default_output_filename(".bw"))

  unmapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: $(default_output_filename(".bw"))
      outputEval: |
        ${
          var unmappedRegex = /align\:.*/;
          return parseInt(self[0].contents.match(unmappedRegex)[0].split(" ")[1]);
        }

  mapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: $(default_output_filename(".bw"))
      outputEval: |
        ${
          var mappedRegex = /alignment\:.*/;
          return parseInt(self[0].contents.match(mappedRegex)[0].split(" ")[1]);
        }

  total_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: $(default_output_filename(".bw"))
      outputEval: |
        ${
          var totalRegex = /processed\:.*/;
          return parseInt(self[0].contents.match(totalRegex)[0].split(" ")[1]);
        }

baseCommand:
  - bowtie

arguments:
  - valueFrom: |
      ${
        if (inputs.upstream_filelist && inputs.downstream_filelist){
          return "-1";
        }
        return null;
      }
    position: 82
  - valueFrom: |
      ${
        if (inputs.upstream_filelist && inputs.downstream_filelist){
          return "-2";
        }
        return null;
      }
    position: 84
  - valueFrom: |
      ${
        return ' 2> ' + default_output_filename(".bw");
      }
    position: 100000
    shellQuote: false

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bowtie_align"
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
  Tool maps input raw reads files to reference genome using Bowtie.

  `default_output_filename` function returns default name for SAM output and log files. In case when `sam` and
  `output_filename` inputs are not set, default filename will have `.sam` extension but format may not correspond SAM
  specification. To set output filename manually use `output_filename` input. Default output filename is based on
  `output_filename` or basename of `upstream_filelist`, `downstream_filelist` or `crossbow_filelist` file (if array,
  the first file in array is taken). If function is called without argenments and `output_filename` input is set, it
  will be returned from the function.

  For single-end input data any of the `upstream_filelist` or `downstream_filelist` inputs can be used.

  Log filename (`log_file` output) is generated by `default_output_filename` function with ex='.bw'

  `indices_folder` defines folder to contain Bowtie indices. Based on the first found file with `rev.1.ebwt` or
  `rev.1.ebwtl` extension, bowtie index prefix is returned from input's `valueFrom` field.