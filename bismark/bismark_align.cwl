#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2

hints:
  SoftwareRequirement:
    packages:
      bismark:
        specs: [ "http://identifiers.org/biotools/bismark" ]
        version: [ "0.0.2" ]

inputs:

  indices_folder:
    type: Directory
    inputBinding:
      position: 3
    label: "Bismark indices folder"
    doc: "Path to Bismark generated indices folder"

  fastq_file:
    type: File
    inputBinding:
      position: 4
    label: "FASTQ file"
    doc: "Uncompressed or gzipped FASTQ file, single-end"

  processes:
    type: int?
    inputBinding:
      position: 1
      prefix: "--multicore"
    label: "Number of Bismark instances to run"
    doc: "Set the number of parallel Bismark instances to run concurrently. Each Bismark instance runs four Bowtie2 aligners"

  threads:
    type: int?
    inputBinding:
      position: 2
      prefix: "-p"
    label: "Number of Bowtie2 threads to use"
    doc: "Set the number of threads for each Bowtie2 aligner"

outputs:

  bam_file:
    type: File
    label: "BAM alignment file"
    doc: "Bismark generated BAM alignment file"
    outputBinding:
      glob: "*.bam"

  alignment_report:
    type: File
    label: "Bismark alignment and methylation report"
    doc: "Bismark generated alignment and methylation summary report"
    outputBinding:
      glob: "*.txt"

baseCommand: ["bismark", "--non_directional"]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "bismark_align"
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
  Default aligner - Bowtie2.
  Only Single-End supported.
  Parameters used:
  --non_directional
    The sequencing library was constructed in a non strand-specific manner, alignments to all four bisulfite strands
    will be reported. (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary
    to the original strands are merely theoretical and should not exist in reality. Specifying directional alignments
    (which is the default) will only run 2 alignment threads to the original top (OT) or bottom (OB) strands in parallel
    and report these alignments. This is the recommended option for strand-specific libraries).