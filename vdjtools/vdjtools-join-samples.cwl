#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool

hints:
  ResourceRequirement:
    ramMin: 3814
    coresMin: 2
  DockerRequirement:
    dockerPull: yyasumizu/vdjtools
requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - var get_label = function(i) {
          var rootname = inputs.molecule_info_h5[i].basename.split('.').slice(0,-1).join('.');
          rootname = (rootname=="")?inputs.molecule_info_h5[i].basename:rootname;
          return inputs.gem_well_labels?inputs.gem_well_labels[i].replace(/,/g, "_"):rootname;
      };
  InitialWorkDirRequirement:
    listing: |
      ${
        var entry = "file.name\tsample.id\n"
        for (var i=0; i < inputs.vdj_file.length; i++){
          entry += inputs.vdj_file[i].path + "\t" + inputs.vdj_name[i] + "\n"
        }
        return [{
          "entry": entry,
          "entryname": "metadata.tsv"
        }];
      }

inputs:
  vdj_file:
    type: File[]
    doc: |
      VDJ formatted files from the vdjtools-convert-mixcr.cwl tool

  vdj_name:
    type: string[]
    doc: |
      Unique names for the files provided in vdj_file

  intersect_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "strict"
      - "nt"
      - "ntV"
      - "ntVJ"
      - "aa"
      - "aaV"
      - "aaVJ"
      - "aa!nt"
    inputBinding:
      prefix: "--intersect-type"
      position: 5
    doc: |
      Specifies which clonotype features (CDR3 sequence, V/J segments, hypermutations)
      will be compared when checking if two clonotypes match.
      Default: aa

  min_times_detected:
    type: int?
    inputBinding:
      prefix: "--times-detected"
      position: 6
    doc: |
      Minimal number of samples where clonotype should be detected.
      Default: 2

  output_prefix:
    type: string?
    default: "./"
    inputBinding:
      position: 7


outputs:
  combined_vdj_file:
    type: File
    outputBinding:
      glob: "*.txt.gz"

  summary_file:
    type: File
    outputBinding:
      glob: "*.summary.txt"

  metadata_file:
    type: File
    outputBinding:
      glob: "metadata.tsv"

  venn_diag_plot:
    type: File
    outputBinding:
      glob: "*.venn.pdf"


baseCommand: ["vdjtools", "JoinSamples", "--compress", "--plot", "--plot-type", "png", "--metadata", "metadata.tsv"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

label: "VDJtools Join Samples"
s:alternateName: "Joins several clonotype tables together to form a joint clonotype abundance table"

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
  VDJtools is an open-source Java/Groovy-based framework designed to facilitate analysis of
  immune repertoire sequencing (RepSeq) data. VDJtools computes a wide set of statistics and
  is able to perform various forms of cross-sample analysis. Both comprehensive tabular output
  and publication-ready plots are provided.
  
  Joins several clonotype tables together to form a joint clonotype abundance table. Joint
  clonotype holds information on all clonotypes that match under a certain comparison criteria
  (e.g. identical CDR3nt and V segment), their samples of origin and corresponding abundances.
