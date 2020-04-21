#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: biocontainers/bandage:v0.8.1-1-deb_cv1

requirements:
  EnvVarRequirement:
    envDef:
      XDG_RUNTIME_DIR: $(runtime.tmpdir)
      QT_QPA_PLATFORM: minimal

label: Bandage image
doc: |
  an hybrid assembly pipeline for bacterial genomes
  *Bandage Overview**
  Bandage is a GUI program that allows users to interact with the assembly graphs made by de novo assemblers 
  such as Velvet, SPAdes,   MEGAHIT and others.
  De novo assembly graphs contain not only assembled contigs but also the connections between those contigs, 
  which were previously not easily accessible. Bandage visualises assembly graphs, with connections, using graph layout algorithms. 
  Nodes in the drawn graph, which represent contigs, can be automatically labelled with their ID, length or depth. Users can interact 
  with the graph by moving, labelling and colouring nodes. Sequence information can also be extracted directly from the graph viewer. 
  By displaying connections between contigs, Bandage opens up new possibilities for analysing and improving de novo assemblies 
  that are not possible by looking at contigs alone. 
  Bandage works with Graphical Fragment Assembly (GFA) files. 
  For more information about this file format, see https://gfa-spec.github.io/GFA-spec/GFA2.html

baseCommand: [ Bandage, image ]

inputs:
  graph:
    type: File
    doc: |
        Graphical Fragment Assembly
        Supports multiple assembly graph formats: 
        LastGraph (Velvet), FASTG (SPAdes), Trinity.fasta, ASQG and GFA.
    inputBinding:
      position: 1

  format:
    type:
      - 'null'
      - type: enum
        symbols:
          - jpg
          - png
          - svg
    default: jpg
    inputBinding:
      position: 2
      valueFrom: $(inputs.graph.nameroot).$(self) 
    doc: |
        Produce jpg, png or svg file

  height:
    type: int
    default: 1000
    inputBinding:
      prefix: --height
      position: 3
    doc: |
        Image height.If only height or width is set,
        the other will be determined automatically.
        If both are set, the image will be exactly that size.

  width:
    inputBinding:
      prefix: --width
      position: 3
    type:  int?
    doc: |
        Image width. If only height or width is set, the other will be determined automatically.
        If both are set, the image will be exactly that size.

  node_length:
    type:  boolean?
    default: true
    inputBinding:
      prefix: --names
      valueFrom: --lengths
      position: 3
    doc: |
        If true, define Node labels as lengths

outputs:
 image:
   type: File
   outputBinding:
      glob: $(inputs.graph.nameroot).$(inputs.format)
   doc: "Assembly Graph Image"
