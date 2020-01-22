#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2

inputs:

  peak_file:
    type: File
    doc: "Homer generated peak file or BED"

  tag_folders:
    type:
      - Directory
      - Directory[]
    inputBinding:
      position: 7
      prefix: "-d"
    doc: "Tag folders from homer-make-tag-directory tool"

  hist_width:
    type:
      - int
      - string
    inputBinding:
      position: 8
      prefix: "-size"
    doc: |
      Possible values:
        "#" - performs analysis on the "#" bp surrounding the peak centers
        "#,#" - performs analysis from "#" to "#" relative to peak center
        "given" - set size to actual coordinates in peak/BED file

  hist_bin_size:
    type: int
    inputBinding:
      position: 9
      prefix: "-hist"
    doc: |
      Bin size, bp. If hist_width is "given" or skipped, this
      parameter will set the number of bins to divide each region into

  export_heatmap:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-ghist"
    doc: |
      Generate heatmap. Instead of averaging all of the data
      from each peak, keep data from each peak separate

  norm_fpkm:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-fpkm"
    doc: |
      Normalize read counts to million reads or fragments per kilobase mapped

  norm_raw:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "-raw"
    doc: |
      Do not adjust the tag counts based on total tags sequenced.
      By default all tag counts will be normalized to norm_tag_count

  norm_tag_count:
    type: int?
    inputBinding:
      position: 13
      prefix: "-norm"
    doc: |
      Normalize tags to this tag count, default=1e7, 0=average tag count in all directories

  norm_fragment_size:
    type: int?
    inputBinding:
      position: 14
      prefix: "-normLength"
    doc: |
      Fragment length to normlize to for experiments with different lens. Default: 100bp

  strand:
    type: string?
    inputBinding:
      position: 15
      prefix: "-strand"
    doc: |
      Count tags on specific strands relative to peak. Default: both
      Possible values: +|-

  threads:
    type: int?
    inputBinding:
      position: 16
      prefix: "-cpu"
    doc: |
      Set the number of threads. Default: 1

  histogram_filename:
    type: string
    doc: "Output histogram's filename"

outputs:

  histogram_file:
    type: stdout
    doc: "Output histogram file"

stdout: ${return inputs.histogram_filename;}

baseCommand: ["annotatePeaks.pl"]
arguments:
  - valueFrom: $(inputs.peak_file)
    position: 5
  - valueFrom: $("none")
    position: 6

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/latest/schema.rdf

s:name: "homer-annotate-peaks-hist"
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
  Tool is used to produce histogram or heatmaps only. Rest of the functionality is not implemented intentionally.
  If TSS analysis needed, input peak_file should be centered on TSS, where the 'center' of the peak in the actual TSS.
  For example:
    1	chr4	978796	978796	-
    2	chr4	1052109	1052109	+
    3	chr4	1105422	1105422	-

  Skipped arguments:

    Related to peaks annotation:
      -organism
      -gtf
      -gff
      -gff3
      -gid
      -ann
      -mask
      -p
      -pdist
      -pcount
      -vcf
      -editDistance
      -individuals
      -gene
      -go
      -genomeOntology
      -ratio
      -rlog
      -vst
      -CpG
      -nfr
      -nfrSize
      -gwasCatalog
      -map
      -noann

    Related to tss/tts/rna modes:
      tss
      tts
      rna
      -list
      -cTSS

    Related to motifs:
      -m
      -mscore
      -nmotifs
      -mdist
      -mfasta
      -fm
      -rmrevopp
      -matrix
      -mbed
      -mlogic
      -norevopp

    Related to peak centering:
      -center
      -mirror
      -multi

    Related to genome comparisons
      -cmpGenome
      -cmpLiftover

    Currently not needed functionality:
      -bedGraph
      -wig
      -nuc
      -di
      -histNorm
      -rm
      -log
      -sqrt
      -len
      -pc
      -noblanks
      -homer1
      -homer2