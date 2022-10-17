#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement:
   expressionLib:
    - |
      var get_output_folder_name = function() {
        if (inputs.output_folder_name == ""){
          var root = inputs.genome_fasta_file.basename.split('.').slice(0,-1).join('.');
          return (root == "") ? inputs.genome_fasta_file.basename : root;
        } else {
          return inputs.output_folder_name;
        }          
      };
  InitialWorkDirRequirement:
    listing: |
      ${
        var exclude_chr = "[]";
        if (inputs.exclude_chr && inputs.exclude_chr.length > 0){
          exclude_chr = '["' + inputs.exclude_chr.join('", "') + '"]'
        }
        var entry = `
        {
            genome: ["${get_output_folder_name()}"]
            input_fasta: ["${inputs.genome_fasta_file.path}"]
            input_gtf: ["${inputs.annotation_gtf_file.path}"]
            non_nuclear_contigs: ${exclude_chr}
        }`
        return [{
          "entry": entry,
          "entryname": "config.txt"
        }];
      }

hints:
  DockerRequirement:
    dockerPull: cumulusprod/cellranger-arc:2.0.0

inputs:
  
  genome_fasta_file:
    type: File
    format: edam:format_1929  # FASTA
    doc: "Genome FASTA file"

  annotation_gtf_file:
    type: File
    format: edam:format_2306  # GTF
    doc: "GTF annotation file."

  exclude_chr:
    type:
      - "null"
      - string[]
    doc: |
      Contigs that do not have any chromatin structure, for example,
      mitochondria or plastids. These contigs are excluded from peak
      calling since the entire contig will be "open" due to a lack of
      chromatin structure

  output_folder_name:
    type: string?
    default: ""
    doc: |
      Unique genome name, used to name output folder

  threads:
    type: int?
    inputBinding:
      position: 5
      prefix: "--nthreads"
    doc: |
      Number of threads used during STAR genome indexing
      Default: 1

  memory_limit:
    type: int?
    inputBinding:
      position: 6
      prefix: "--memgb"
    doc: |
      Maximum memory (GB) used when aligning reads with STAR
      Defaults: 16


outputs:
  indices_folder:
    type: Directory
    outputBinding:
      glob: $(get_output_folder_name())
    doc: |
      Compatible with Cell Ranger ARC reference folder that includes
      STAR and BWA indices

baseCommand: ["cellranger-arc", "mkref", "--config", "config.txt"]

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger ARC mkref - builds compatible with Cell Ranger ARC indices"
s:alternateName: "Builds compatible with Cell Ranger ARC reference folder from user-supplied genome FASTA and gene GTF files"

s:license: http://www.apache.org/licenses/LICENSE-2.0

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
  Reference preparation tool for 10x Genomics Cell Ranger Multiome ATAC + Gene Expression

  Notes:
  - `input_motifs` parameter in the `config.txt` file is not implemented
  - if GTF file provided in `annotation_gtf_file` has duplicate gene_id, they should be
    grouped together. Applicable to to USCS RefGene annotations.
