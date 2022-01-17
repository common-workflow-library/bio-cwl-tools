cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.fasta_file.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.fasta_file.basename+".2bit":root+".2bit";
          } else {
            return inputs.output_filename;
          }
        };


inputs:

  fasta_file:
    type: File
    inputBinding:
      position: 5
    doc: "Reference genome FASTA file"

  output_filename:
    type: string?
    default: ""
    inputBinding:
      valueFrom: $(default_output_filename())
      position: 6
    doc: "Output file name"


outputs:

  twobit_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Reference genome 2bit file"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["faToTwoBit"]
stdout: fa_to_twobit_stdout.log
stderr: fa_to_twobit_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "ucsc-fa-to-twobit"
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
  faToTwoBit - Convert DNA from fasta to 2bit format