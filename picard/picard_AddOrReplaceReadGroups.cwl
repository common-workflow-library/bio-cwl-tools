#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand:
- picard
- AddOrReplaceReadGroups

doc: |-
  Assigns all the reads in a file to a single new read-group.

   <h3>Summary</h3>
   Many tools (Picard and GATK for example) require or assume the presence of at least one <code>RG</code> tag, defining a "read-group"
   to which each read can be assigned (as specified in the <code>RG</code> tag in the SAM record).
   This tool enables the user to assign all the reads in the INPUT to a single new read-group.
   For more information about read-groups, see the <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>
   GATK Dictionary entry.</a>
   <br />
   This tool accepts as INPUT BAM and SAM files or URLs from the
   <a href="http://ga4gh.org/#/documentation">Global Alliance for Genomics and Health (GA4GH)</a>.
   <h3>Caveats</h3>
   The value of the tags must adhere (according to the <a href="https://samtools.github.io/hts-specs/SAMv1.pdf">SAM-spec</a>)
   with the regex <pre>#READGROUP_ID_REGEX</pre> (one or more characters from the ASCII range 32 through 126). In
   particular <code>&lt;Space&gt;</code> is the only non-printing character allowed.
   <br/>
   The program enables only the wholesale assignment of all the reads in the INPUT to a single read-group. If your file
   already has reads assigned to multiple read-groups, the original <code>RG</code> value will be lost.
  Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement:
    expressionLib:
    - |
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
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/picard:2.22.2--0
inputs:
- doc: Input file (BAM or SAM or a GA4GH url). [synonymous with -I]
  id: INPUT
  type: File
  inputBinding:
    prefix: INPUT=
    separate: false
- doc: Read-Group library [synonymous with -LB]
  id: RGLB
  type: string
  inputBinding:
    prefix: RGLB=
    separate: false
- doc: Read-Group platform (e.g. ILLUMINA, SOLID) [synonymous with -PL]
  id: RGPL
  type: string
  inputBinding:
    prefix: RGPL=
    separate: false
- doc: Read-Group platform unit (eg. run barcode) [synonymous with -PU]
  id: RGPU
  type: string
  inputBinding:
    prefix: RGPU=
    separate: false
- doc: Read-Group sample name [synonymous with -SM]
  id: RGSM
  type: string
  inputBinding:
    prefix: RGSM=
    separate: false
- doc: Output filename (BAM or SAM)
  id: OUTPUT
  type: string
  inputBinding:
    prefix: OUTPUT=
    separate: false
- doc: Reference sequence file. [synonymous with -R]
  id: REFERENCE_SEQUENCE
  type: File?
  inputBinding:
    prefix: REFERENCE_SEQUENCE=
    separate: false
- doc: Optional sort order to output in. If not supplied OUTPUT is in the same order
    as INPUT. [synonymous with -SO]
  id: SORT_ORDER
  type:
  - 'null'
  - type: enum
    symbols:
    - unsorted
    - queryname
    - coordinate
    - duplicate
    - unknown
  inputBinding:
    prefix: SORT_ORDER=
    separate: false
- doc: Read-Group sequencing center name [synonymous with -CN]
  id: RGCN
  type: string?
  inputBinding:
    prefix: RGCN=
    separate: false
- doc: Read-Group description [synonymous with -DS]
  id: RGDS
  type: string?
  inputBinding:
    prefix: RGDS=
    separate: false
- doc: Read-Group run date in Iso8601Date format [synonymous with -DT]
  id: RGDT
  type: string?
  inputBinding:
    prefix: RGDT=
    separate: false
- doc: Read-Group flow order [synonymous with -FO]
  id: RGFO
  type: string?
  inputBinding:
    prefix: RGFO=
    separate: false
- doc: Read-Group ID [synonymous with -ID]
  id: RGID
  type: string?
  inputBinding:
    prefix: RGID=
    separate: false
- doc: Read-Group key sequence [synonymous with -KS]
  id: RGKS
  type: string?
  inputBinding:
    prefix: RGKS=
    separate: false
- doc: Read-Group program group [synonymous with -PG]
  id: RGPG
  type: string?
  inputBinding:
    prefix: RGPG=
    separate: false
- doc: Read-Group predicted insert size [synonymous with -PI]
  id: RGPI
  type: int?
  inputBinding:
    prefix: RGPI=
    separate: false
- doc: Read-Group platform model [synonymous with -PM]
  id: RGPM
  type: string?
  inputBinding:
    prefix: RGPM=
    separate: false
- doc: Control verbosity of logging.
  id: VERBOSITY
  type:
  - 'null'
  - type: enum
    symbols:
    - ERROR
    - WARNING
    - INFO
    - DEBUG
  inputBinding:
    prefix: VERBOSITY=
    separate: false
- doc: Whether to suppress job-summary info on System.err.
  id: QUIET
  type: boolean?
  inputBinding:
    prefix: QUIET=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: Validation stringency for all SAM files read by this program.  Setting stringency
    to SILENT can improve performance when processing a BAM file in which variable-length
    data (read, qualities, tags) do not otherwise need to be decoded.
  id: VALIDATION_STRINGENCY
  type:
  - 'null'
  - type: enum
    symbols:
    - STRICT
    - LENIENT
    - SILENT
  inputBinding:
    prefix: VALIDATION_STRINGENCY=
    separate: false
- doc: Compression level for all compressed files created (e.g. BAM and VCF).
  id: COMPRESSION_LEVEL
  type: int?
  inputBinding:
    prefix: COMPRESSION_LEVEL=
    separate: false
- doc: When writing files that need to be sorted, this will specify the number of
    records stored in RAM before spilling to disk. Increasing this number reduces
    the number of file handles needed to sort the file, and increases the amount of
    RAM needed.
  id: MAX_RECORDS_IN_RAM
  type: int?
  inputBinding:
    prefix: MAX_RECORDS_IN_RAM=
    separate: false
- doc: Use the JDK Deflater instead of the Intel Deflater for writing compressed output
    [synonymous with -use_jdk_deflater]
  id: USE_JDK_DEFLATER
  type: boolean?
  inputBinding:
    prefix: USE_JDK_DEFLATER=
    separate: false
    valueFrom: $(generateGATK4BooleanValue())
- doc: Use the JDK Inflater instead of the Intel Inflater for reading compressed input
    [synonymous with -use_jdk_inflater]
  id: USE_JDK_INFLATER
  type: boolean?
  inputBinding:
    prefix: USE_JDK_INFLATER=
    separate: false
    valueFrom: $(generateGATK4BooleanValue())
- doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
  id: CREATE_INDEX
  type: boolean?
  inputBinding:
    prefix: CREATE_INDEX=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: 'Whether to create an MD5 digest for any BAM or FASTQ files created.  '
  id: CREATE_MD5_FILE
  type: boolean?
  inputBinding:
    prefix: CREATE_MD5_FILE=
    valueFrom: $(generateGATK4BooleanValue())
    separate: false
- doc: Google Genomics API client_secrets.json file path.
  id: GA4GH_CLIENT_SECRETS
  type: File?
  inputBinding:
    prefix: GA4GH_CLIENT_SECRETS=
    separate: false

arguments:
 - TMP_DIR=$(runtime.tmpdir)

outputs:
  sequences_with_new_read_group:
    type: File
    format: edam:format_2573  # SAM
    outputBinding:
      glob: $(inputs.OUTPUT)

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
