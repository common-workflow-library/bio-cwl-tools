cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            let root = inputs.bam_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.bam_file.basename+".bed":root+".bed";
          } else {
            return inputs.output_filename;
          }
        };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0


inputs:

  bam_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "Input BAM file (not necessary sorted or indexed)"

  output_filename:
    type: string?
    default: ""
    doc: "Output BED filename"


outputs:

  bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Sequences file"


baseCommand: ["bedtools", "bamtobed"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/bedtools_metadata.yaml

s:name: "bedtools_bamtobed_2"
s:downloadUrl: https://raw.githubusercontent.com/common-workflow-library/bio-cwl-tools/release/tools/bedtools/bedtools_bamtobed_2.cwl
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
  Converts BAM to BED. All Options are not implemented.

s:about: |
  Usage:   bedtools bamtobed [OPTIONS] -i <bam> 

  Options: 
    -bedpe	Write BEDPE format.
      - Requires BAM to be grouped or sorted by query.

    -mate1	When writing BEDPE (-bedpe) format, 
      always report mate one as the first BEDPE "block".

    -bed12	Write "blocked" BED format (aka "BED12"). Forces -split.

      http://genome-test.cse.ucsc.edu/FAQ/FAQformat#format1

    -split	Report "split" BAM alignments as separate BED entries.
      Splits only on N CIGAR operations.

    -splitD	Split alignments based on N and D CIGAR operators.
      Forces -split.

    -ed	Use BAM edit distance (NM tag) for BED score.
      - Default for BED is to use mapping quality.
      - Default for BEDPE is to use the minimum of
        the two mapping qualities for the pair.
      - When -ed is used with -bedpe, the total edit
        distance from the two mates is reported.

    -tag	Use other NUMERIC BAM alignment tag for BED score.
      - Default for BED is to use mapping quality.
        Disallowed with BEDPE output.

    -color	An R,G,B string for the color used with BED12 format.
      Default is (255,0,0).

    -cigar	Add the CIGAR string to the BED entry as a 7th column.