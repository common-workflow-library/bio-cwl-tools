cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.genome_folder,
                  "writable": true
                }
              ]
    }

inputs:

  genome_folder:
    type: Directory
    inputBinding:
      position: 2
    label: "Genome folder"
    doc: "Genome folder with FASTA files"


outputs:

  indices_folder:
    type: Directory
    label: "Bismark indices folder"
    doc: "Bismark generated indices folder"
    outputBinding:
      glob: "*"


baseCommand: ["bismark_genome_preparation"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bismark_metadata.yaml

s:name: "bismark_prepare_genome"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bismark/bismark_prepare_genome.cwl
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
  bismark_genome_preparation script generates indices using Bowtie2 aligners by default.


s:about: |
  USAGE: bismark_genome_preparation [options] <arguments>

  OPTIONS:

  --help/--man             Displays this help filea and exits.

  --version                Displays version information and exits.

  --verbose                Print verbose output for more details or debugging.

  --path_to_bowtie </../>  The full path to the Bowtie 1 or Bowtie 2 installation on your system
                          (depending on which aligner/indexer you intend to use). Unless this path
                          is specified it is assumed that Bowtie is in the PATH.

  --bowtie2                This will create bisulfite indexes for Bowtie 2. (Default: ON).

  --bowtie1                This will create bisulfite indexes for Bowtie 1. (Default: OFF).

  --single_fasta           Instruct the Bismark Indexer to write the converted genomes into
                          single-entry FastA files instead of making one multi-FastA file (MFA)
                          per chromosome. This might be useful if individual bisulfite converted
                          chromosomes are needed (e.g. for debugging), however it can cause a
                          problem with indexing if the number of chromosomes is vast (this is likely
                          to be in the range of several thousand files; the operating system can
                          only handle lists up to a certain length, and some newly assembled
                          genomes may contain 20000-500000 contigs of scaffold files which do exceed
                          this list length limit).

  --genomic_composition    Calculate and extract the genomic sequence composition for mono and di-nucleotides
                          and write the genomic composition table 'genomic_nucleotide_frequencies.txt' to the
                          genome folder. This may be useful later on when using bam2nuc or the Bismark option
                          --nucleotide_coverage.

  ARGUMENTS:

  <path_to_genome_folder>  The path to the folder containing the genome to be bisulfite converted.
                          The Bismark Genome Preparation expects one or more fastA files in the folder
                          (with the file extension: .fa or .fasta). Specifying this path is mandatory.