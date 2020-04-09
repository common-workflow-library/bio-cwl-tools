cwlVersion: v1.0
class: CommandLineTool
id: bandage-info
inputs:

  - id: graph
    type:  File
    doc: |
        Graphical Fragment Assembly.
        Supports multiple
        assembly graph formats: 
        LastGraph (Velvet), FASTG (SPAdes), Trinity.fasta, ASQG and GFA.    


  - id: tsv
    type:  boolean
    default: false
    doc: |
        If true, output the information in a single tab-delimited line 
        starting with the graph file





outputs:

 - id: all_script
   type:
      - type: array
        items: File
   outputBinding:
      glob: "*.sh"  
   doc: "generated script to run bandage. for learning purpose" 

 - id: assembly_graph_info
   type: File
   outputBinding:
      glob: "assembly_graph_info.txt"
   doc: "Assembly Graph Information"




baseCommand: bash

arguments: [bandage_info_launch.sh]

hints:
  DockerRequirement:
    dockerPull: "fjrmore/bandage" 


requirements:
  - class:  InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: bandage_info_launch.sh
        entry: |
               #!/bin/bash
               ###########################
               # Bandage info wrapper  
               export QT_QPA_PLATFORM=offscreen
               TMPDIR=$PWD"/tmp_runtime-bandage"
               mkdir -p $TMPDIR
               export XDG_RUNTIME_DIR=$TMPDIR
               Bandage info '$(inputs.graph.path)' \\
               ${
                var opt=""
                if(inputs.tsv==true){ 
                 opt+=" --tsv "
                }
                return opt
               } \\
                > assembly_graph_info.txt


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

