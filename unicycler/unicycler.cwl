class: CommandLineTool
cwlVersion: v1.0

id: unicycler
baseCommand:
  - bash
#  - cat
inputs:


 
##tyoe selector
  - id: fastq_file_type
    type:
      type: enum
      symbols: [paired,single]
#     symbols: [paired,single,paired_collection]
#    default: single   
    doc:  "Paired and single end data"

##############

##input1

  - id: fastq1_type
    type:
      type: enum
      symbols:
        - fastqsanger
        - fastqsanger.gz
    doc: "Type of the First set of reads. Only when fastq_file_type = single  or  paired"  
    default: fastqsanger

  - id: fastq1
    type: File
    doc: "First set of reads with forward reads. Only when fastq_file_type = single or paired" 



##input2


  - id: fastq2_type
    type:
     - "null"
     -  type: enum
        symbols:
         - fastqsanger
         - fastqsanger.gz
    default: "null"
    doc: "Type of the Second set of reads. Only when fastq_file_type=paired"  

    
  - id: fastq2
    type: File?
    doc: "Second set of reads with reverse reads. Only when fastq_file_type=paired"  

 





#sequence_long (optional)

  - id: sequence_long_type
    type:
     - "null"
     -  type: enum
        symbols:
          - fastqsanger
          - fastqsanger.gz
          - fasta
    doc: "long reads file type. If there are no long reads, leave this empty"

    
  - id: sequence_long
    type:  File?
    doc: "long reads. If there are no long reads, leave this empty"
#    
  - id: compute_slots
    type:  int
    default: 4
    doc: ""
    



##mode selector

  - id: mode
    type:
      type: enum
      symbols: [conservative,normal,bold]
    doc: |
        Bridging mode, values:
        conservative (smaller contigs, lower misassembly)
        normal (moderate contig size and misassembly rate)
        bold  (longest contigs, higher misassembly rate)


  - id: min_fasta_length
    type:  int
    default: 100
    doc: "Exclude contigs from the FASTA file which are shorter than this length (bp)"


  - id: linear_seqs
    type:  int
    default: 0
    doc: "The expected number of linear (i.e. non-circular) sequences in the assembly"


  - id: min_anchor_seg_len
    type:  int
    default: 0
    doc: "Unicycler will not use segments shorter than this as scaffolding anchors"


  - id: spades_no_correct
    type:  boolean
    default: false
    doc: |
        Unicycler uses SPAdes to construct assembly graphs.
        You can modify some of the SPAdes settings here.
        Use this ONLY if you know what you are doing!
        This option turns off SPAdes error correction. 
        Generally it is highly recommended to use correction.

 
  - id: spades_min_kmer_frac
    type:  float
    default: 0.2
    doc: |
         Lowest k-mer size for SPAdes assembly,
         expressed as a fraction of the read length. 
         min 0, max 1

  - id: spades_max_kmer_frac
    type:  float
    default: 0.95
#    max: 1.0
#    min: 0.0
    doc: |
         Highest k-mer size for SPAdes assembly,
         expressed as a fraction of the read length. 
         min 0, max 1

  
  # check default values
  - id: spades_kmers
    type:  string
    default: "11,127"
    doc: |
         Exact k-mers size to use for SPAdes assembly, comma-separated"
         Kmers must be comma-separated odd integers (no repitition)
         without space in the range of 11 to 127 (inclusive)

  - id: spades_kmer_count
    type:  int
    default: 10
#    min: 0
    doc: "Number of k-mer steps to use in SPAdes assembly, min 0"


 
  - id: spades_depth_filter
    type:  float
    default: 0.25
#    min: 0.0
#    max: 1.0
    doc: |
         Filter out contigs lower than this fraction 
         of the chromosomal depth.
         It is done if does not result in graph dead ends
         min 0, max 1

  - id: spades_largest_component
    type:  boolean
    default: false
    doc: "Only keep the largest connected component of the assembly graph if true"



  - id: rotation_no_rotate
    type:  boolean
    default: false
    doc: |
         These options control the rotation of completed circular sequence
         near the end of the Unicycler pipeline. Use this ONLY if you know what you are doing!
         Do not rotate completed replicons to start at a standard gene.
         Unicycler uses TBLASTN to search for dnaA or repA alleles in each completed replicon.
         If one is found, the sequence is rotated and/or flipped so that it begins with that gene
         encoded on the forward strand. This provides consistently oriented assemblies and reduces
         the risk that a gene will be split across the start and end of the sequence.
    
  - id: rotation_start_genes
    type:  File? 
    doc: "FASTA file of genes for start point of rotated replicons"



  - id: rotation_start_gene_id
    type:  float
    default: 90.0
    doc: "The minimum required BLAST percent identity for a start gene search. max 100, min 0"


  - id: rotation_start_gene_cov
    type:  float
    default: 95.0
    doc: "The minimum required BLAST percent coverage for a start gene search. min 0, max 100"


  - id: pilon_no_pilon
    type:  boolean
    default: false
    doc: "Unicycler uses Pilon tool for polishing final assembly. Do not use if true"


  - id: graph_clean_min_component_size
    type:  int
    default: 1000
    doc: "Contigs shorter than this value (bp) will not be polished using Pilon; min 0"

    
  - id: graph_clean_min_dead_end_size
    type:  int
    default: 1000
    doc: |
         These options control the removal of small leftover sequences after bridging is complete.
         Unbridged graph components smaller than this size will be removed from the final graph,
         min 0    

  - id: lr_align_contamination
    type:  File?
    doc: |
         FASTA file of known contamination in long reads,
         e.g. lambda, phiXm or puc18 spike-ins. 

  - id: lr_align_scores
    type:  string
    default: "3,-6,-5,-2"
    doc: |
          Comma-delimited string of alignment scores: match, mismatch, gap open, gap extend


  - id: lr_align_low_score
    type:  int?
    doc: |
          Score threshold - alignments below this are considered poor,
          default = set automatically



 
    
outputs:
  - id: exec_script
    type: File
    outputBinding:
      glob: "*.sh"  
    doc: |
          Launching script for learning purpose


  - id: assembly_graph
    type: File
    outputBinding:
      glob: assembly.gfa
    doc: |
          Assembly graph



  - id: assembly
    type: File
    outputBinding:
      glob: assembly.fasta
    doc: |
          fasta assembly output sequence
          (main output)
      
      
arguments:
  - unicycler_launch.sh
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: unicycler_launch.sh
        entry: |
               #!/bin/bash
               ###########################
               #      unicycler launcher
               ###########################

               ##preparing input files
               #check permission / chmod  is issues
               ${
                 var fl=""
                 var lncmd="" 
                 var fq1=""
                 var fq2=""
                 var lr=""
                  
               //###################paired case              
                     if (inputs.fastq_file_type =="paired"  ){
                      if( inputs.fastq1_type=='fastqsanger' ){
                          fq1 = "fq1.fastq"
                      }
                      else if( inputs.fastq1_type=='fastqsanger.gz' ){
                           fq1 = "fq1.fastq.gz"
                      }
                      if( inputs.fastq2_type=='fastqsanger' ){
                          fq2 = "fq2.fastq"
                       }
                       else if( inputs.fastq2_type=='fastqsanger.gz' ){
                           fq2 = "fq2.fastq.gz"
                        }
                        lncmd+="fq1='"+fq1+"'"
                        lncmd+=" && "
                        lncmd+="fq2='"+fq2+"'"
                        lncmd+=" && "
                        lncmd+=" ln -s '"+inputs.fastq1.path+"' $fq1 "
                        lncmd+=" && "
                        lncmd+=" ln -s '"+inputs.fastq2.path+"' $fq2  "
        
                     }
                //###################single case
                 
                if (inputs.fastq_file_type =="single"  ){
                  if( inputs.fastq1_type=='fastqsanger' ){
                      fq1 = "fq1.fastq"
                  }
                  else if( inputs.fastq1_type=='fastqsanger.gz' ){
                       fq1 = "fq1.fastq.gz"
                  }
                  lncmd+="fq1='"+fq1+"'"
                  lncmd+=" && "
                  lncmd+=" ln -s '"+inputs.fastq1.path+"' $fq1 "
                 }
                 //####### long reads
                  if (  inputs.sequence_long !== null) {
                      if (inputs.sequence_long_type=='fastqsanger'){
                               lr = "lr.fastq"
                      }
                      else if (inputs.sequence_long_type=='fastqsanger.gz') {
                               lr = "lr.fastq.gz"
                      }
                      else if (inputs.sequence_longg_type=='fasta') {
                               lr = "lr.fasta"
                      }
                      lncmd+="lr='"+lr+"'"
                      lncmd+=" && "
                      lncmd+= " ln -s '"+inputs.sequence_long.path+"' '$lr' "
                  }


                  return lncmd

               }


               ##general options

               read -d '' GENERALOPT << EOF
               ${
                var opt=""  
                //## General Unicycler Options section
                opt+=" --mode "+inputs.mode+" "
                opt+=" --min_fasta_length "+inputs.min_fasta_length+" "
                opt+=" --linear_seqs "+inputs.linear_seqs+" "

                if (inputs.min_anchor_seg_len  != null ){opt+=" --min_anchor_seg_len "+inputs.min_anchor_seg_len+" "}

                //## Spades Options section
                if(inputs.spades_no_correct==true){opt+=" --no_correct "}
                opt+=" --min_kmer_frac "+inputs.spades_min_kmer_frac+" "
                opt+=" --max_kmer_frac "+inputs.spades_max_kmer_frac+" "
                if (inputs.spades_kmers   != null){opt+=" --kmers "+inputs.spades_kmers+" "}
                               
                opt+=" --kmer_count "+inputs.spades_kmer_count+" "
                opt+=" --depth_filter "+inputs.spades_depth_filter+" "
                if (inputs.spades_largest_component){opt+=" --largest_component "}
                //## Rotation Options section
                if(inputs.rotation_no_rotate == true){ opt+=" --no_rotate "}
                if (inputs.rotation_start_genes!=null){opt+=" --start_genes "+ inputs.rotation_start_genes.path+ " "}
                opt+=" --start_gene_id "+inputs.rotation_start_gene_id+" "
                opt+=" --start_gene_cov "+inputs.rotation_start_gene_cov+" "
                return opt
                }
               EOF
               
               ##additionnal option
              
               read -d '' ADDOPT << EOF
               ${

                var opt=""
                
                if (inputs.pilon_no_pilon  == true){ opt+=" --no_pilon " }
                if (inputs.pilon_min_polish_size  != null){opt+=" --min_polish_size "+inputs.pilon_min_polish_size + " "}
                //## Long Read Alignment Options
                if ( inputs.lr_align_contamination!=null){opt+=" --contamination "+inputs.lr_align_contamination + " "}
                opt+=" --scores "+inputs.lr_align_scores+" "
                if (inputs.lr_align_low_score != null){opt+=" --low_score "+inputs.lr_align_low_score+" "}
                 return ''+ opt + ''
               }
               EOF

               ## Get location for pilon jar file

               ${
                 var cmd=""  
                 cmd+="PILONJAR=/usr/share/java/pilon.jar "
                 return cmd
               }   

               ## Build Unicycler command
               ${
               
                 var cmd_base=""
                 var opt=""

                 
                 
                 
                 cmd_base+=" unicycler -t "+inputs.compute_slots+"  "
                 cmd_base+=" -o ./  "
                 cmd_base+=" --verbosity 3  "
                 cmd_base+=" --pilon_path \$PILONJAR  "

                
                if ( inputs.fastq_file_type == "paired"){
                       opt+=" -1 $fq1 -2 $fq2  "
                }
                else if ( inputs.fastq_file_type == "paired_collection"){
                       opt+=" -1 $fq1 -2 $fq2  "
                }  
                else if ( inputs.fastq_file_type == "single"){
                   opt+=" -s $fq1 "
                }
                if (  inputs.sequence_long !== null) {
                  opt+=" -l $lr "
                }


                //##  Unicycler command
                var cmdl=cmd_base+" "+opt+" \$GENERALOPT \$ADDOPT "
               
                return cmdl

                }
                

 
#        writable: false
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biocontainers/unicycler:v0.4.7dfsg-2-deb_cv1'
  
doc: |
  CWL  wrapped for Unicycler.
  an hybrid assembly pipeline for bacterial genomes
  see  https://github.com/rrwick/Unicycler
  outputs
    final assembly in FASTA format (major output)
    final assembly grapth in graph format, visualized using tools such as  
    Bandage  https://github.com/rrwick/Bandage      


