cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ".fastxstat"
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastx_toolkit:v0.0.14


inputs:

  input_file:
    type: File
    inputBinding:
      position: 10
      prefix: -i
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)

  new_output_format:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '-N'
    doc: |
      New output format (with more information per nucleotide/cycle).
      cycle (previously called 'column') = cycle number
      max-count
      For each nucleotide in the cycle (ALL/A/C/G/T/N):
          count   = number of bases found in this column.
          min     = Lowest quality score value found in this column.
          max     = Highest quality score value found in this column.
          sum     = Sum of quality score values for this column.
          mean    = Mean quality score value for this column.
          Q1	= 1st quartile quality score.
          med	= Median quality score.
          Q3	= 3rd quartile quality score.
          IQR	= Inter-Quartile range (Q3-Q1).
          lW	= 'Left-Whisker' value (for boxplotting).
          rW	= 'Right-Whisker' value (for boxplotting).

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 11
      prefix: -o
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
      Output file to store generated statistics. If not provided - return from default_output_filename function


outputs:

  statistics_file:
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
    doc: Generated statistics file


baseCommand: [fastx_quality_stats]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/fastx_toolkit_metadata.yaml

s:name: "fastx_quality_stats"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/fastx_toolkit/fastx_quality_stats.cwl
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
  Tool calculates statistics on the base of FASTQ file quality scores.
  If `output_filename` is not provided call function `default_output_filename` to return default output file name
  generated as `input_file` basename + `.fastxstat` extension.

s:about: |
  usage: fastx_quality_stats [-h] [-N] [-i INFILE] [-o OUTFILE]
  Part of FASTX Toolkit 0.0.14 by A. Gordon (assafgordon@gmail.com)
     [-h] = This helpful help screen.
     [-i INFILE]  = FASTQ input file. default is STDIN.
     [-o OUTFILE] = TEXT output file. default is STDOUT.
     [-N]         = New output format (with more information per nucleotide/cycle).
  The *OLD* output TEXT file will have the following fields (one row per column):
  	column	= column number (1 to 36 for a 36-cycles read solexa file)
  	count   = number of bases found in this column.
  	min     = Lowest quality score value found in this column.
  	max     = Highest quality score value found in this column.
  	sum     = Sum of quality score values for this column.
  	mean    = Mean quality score value for this column.
  	Q1	= 1st quartile quality score.
  	med	= Median quality score.
  	Q3	= 3rd quartile quality score.
  	IQR	= Inter-Quartile range (Q3-Q1).
  	lW	= 'Left-Whisker' value (for boxplotting).
  	rW	= 'Right-Whisker' value (for boxplotting).
  	A_Count	= Count of 'A' nucleotides found in this column.
  	C_Count	= Count of 'C' nucleotides found in this column.
  	G_Count	= Count of 'G' nucleotides found in this column.
  	T_Count	= Count of 'T' nucleotides found in this column.
  	N_Count = Count of 'N' nucleotides found in this column.
  	max-count = max. number of bases (in all cycles)
  The *NEW* output format:
  	cycle (previously called 'column') = cycle number
  	max-count
    For each nucleotide in the cycle (ALL/A/C/G/T/N):
  		count   = number of bases found in this column.
  		min     = Lowest quality score value found in this column.
  		max     = Highest quality score value found in this column.
  		sum     = Sum of quality score values for this column.
  		mean    = Mean quality score value for this column.
  		Q1	= 1st quartile quality score.
  		med	= Median quality score.
  		Q3	= 3rd quartile quality score.
  		IQR	= Inter-Quartile range (Q3-Q1).
  		lW	= 'Left-Whisker' value (for boxplotting).
  		rW	= 'Right-Whisker' value (for boxplotting).
