cwlVersion: v1.0
class: CommandLineTool
id: bandage-image
inputs:

  - id: graph
    type: File
    doc: |
        Graphical Fragment Assembly
        Supports multiple assembly graph formats: 
        LastGraph (Velvet), FASTG (SPAdes), Trinity.fasta, ASQG and GFA.
 

  - id: format
    type:  string
    default: jpg
    doc: |
        Produce jpg, png or svg file


  - id: height
    type:  int
    default: 1000
    doc: |
        Image height.If only height or width is set, 
        the other will be determined automatically.
        If both are set, the image will be exactly that size.


  - id: width
    type:  int?
    doc: |
        Image width. If only height or width is set, the other will be determined automatically.
        If both are set, the image will be exactly that size.


  - id: node_name
    type:  boolean
    default: true
    doc: |
        If true, define Node labels as name 

  - id: node_length
    type:  boolean
    default: true
    doc: |
        If true, define Node labels as length 


outputs:


 - id: all_script
   type:
      - type: array
        items: File
   outputBinding:
      glob: "*.sh"  
   doc: "generated script to run bandage. for learning purpose" 


 - id: image
   type: File
   outputBinding:
      glob: "*.$(inputs.format)"
   doc: "Assembly Graph Image"




baseCommand: bash

arguments: [bandage_image_launch.sh]

hints:
  DockerRequirement:
    dockerPull: "fjrmore/bandage" 


requirements:
  - class:  InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: bandage_image_launch.sh
        entry: |
               #!/bin/bash
               ###########################
               # Bandage image wrapper  
               export QT_QPA_PLATFORM=minimal
               TMPDIR=$PWD"/tmp_runtime-bandage"
               mkdir -p $TMPDIR
               export XDG_RUNTIME_DIR=$TMPDIR
               GRAPH="$(inputs.graph.path)"
               IMAGE="$(inputs.graph.nameroot).$(inputs.format)"
               Bandage image $GRAPH $IMAGE  \\
               ${
                var opt=""
                if(inputs.height!=null){ 
                 opt+=" --height "+inputs.height+ " "
                }
                if(inputs.width!=null){ 
                 opt+=" --width "+inputs.width +" "
                }
                if(inputs.node_length==true){ 
                 opt+=" --names "
                }
                if(inputs.node_length==true){ 
                 opt+=" --lengths "
                }
                return opt
               }  


doc: |
  CWL  tool for Bandage-info.
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

