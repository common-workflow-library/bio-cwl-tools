#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

baseCommand:
  - bash

inputs:
  libraries_metadata:
    type:
      type: array
      items:
          type: record
          fields:
             lib_index: int? 
             orientation: string?
             lib_type: string? 
    doc: |
        reads library metadata
        related to   libraries_fwd_rev and libraries_mono inputs
        lib_index(id) must match
  libraries_fwd_rev:
    type:
      type: array
      items:
          type: record
          fields:
             lib_index: int? 
             fwd_reads: File?
             rev_reads: File?
    doc: |
        reads file
        orientation must be a value in  ff, fr, rf
        K-mer choices can be chosen by SPAdes instead of being entered manually
  libraries_mono:
    type:
      type: array
      items:
          type: record
          fields:
             lib_index: int? 
             file_type: string?
             reads: File?
    doc: |
        reads file
        file_type value must be in : interleaved, merged, unpaired
  pacbio_reads:
    type: File[]?
  nanopore_reads:
    type: File[]?
  sanger_reads:
    type: File[]?
  trusted_contigs:
    type: File[]?
  untrusted_contigs:
    type: File[]?
  auto_kmer_choice:
    type:  boolean
    default: true
    doc: |
        Automatically choose k-mer values.
        K-mer choices can be chosen by SPAdes instead of being entered manually
  kmers:
    type:  string
    default: "21,33,55"
    doc: |
        K-mers to use, separated by commas.
        Comma-separated list of k-mer sizes to be used 
        (all values must be odd, less than 128, listed in ascending order,
         and smaller than the read length). The default value is 21,33,55
  cov_state:
    type:
     - "null"
     -  type: enum
        symbols:
          - off
          - value
          - auto
    doc: |
        Coverage cutoff ( 'auto', or 'off', or 'value'). auto if null
        when cov_state=value (User Specific) , cov_cutoff must be provided
  cov_cutoff:
    type:  float?
    doc: |
        coverage cutoff value (a positive float number )
  iontorrent:
    type:  boolean
    default: false
    doc: |
        true if Libraries are IonTorrent reads.
  sc:
    type:  boolean
    default: false
    doc: |
        This option is required for MDA. 
        true if single-cell data. 
  onlyassembler:
    type: boolean
    default: false
    doc: |
        Run only assembly if true
        (without read error correction)
  careful:
    type: boolean
    default: true
    doc: |
        Careful correction.
        Tries to reduce number of mismatches and short indels. 
        Also runs MismatchCorrector, a post processing tool,
        which uses BWA tool (comes with SPAdes).

outputs:
  out_contig_stats:
    type: File
    outputBinding:
      glob: out_contig_stats.*
    doc: "contig stats, default column_names: name,length,coverage"
  out_scaffold_stats:
    type: File
    outputBinding:
      glob: out_scaffold_stats.*
    doc: "scaffold stats, default column_names: name,length,coverage"

  out_contigs:
    type: File
    outputBinding:
      glob: contigs.fasta
    doc: "contigs (fasta sequence)"

  out_scaffolds:
    type: File
    outputBinding:
      glob: scaffolds.fasta
    doc: "scaffolds (fasta sequence)"
#  - id: all_script
#    type:
#      - type: array
#        items: File
#    outputBinding:
#      glob: "*.sh"  
#    doc: "generated script to run spades. for learning purpose" 
  all_log:
    type:
      - type: array
        items: File
    outputBinding:
      glob: "*.log"  
    doc: "spades output log and warnings" 
#  - id: out_log
#    type: File
#    outputBinding:
#      glob: spades.log
#    doc: "spades output log"
  assembly_graph:
    type: File
    outputBinding:
      glob: assembly_graph.fastg
    doc: "assembly graph" 
  assembly_graph_with_scaffolds:
    type: File
    outputBinding:
      glob: assembly_graph_with_scaffolds.gfa
    doc: "assembly graph with scaffolds" 
      
arguments:
  - spades_wrapper.sh
requirements:
  - class: InitialWorkDirRequirement
    listing:
####################
      - entryname: write_tsv.py
        entry: |
               #######input file here
               #!/usr/bin/env python
               import sys,re
               search_str = r'^>(NODE|\S+)_(\d+)(?:_|\s)length_(\d+)_cov_(\d+\.*\d*).*\$'
               replace_str = r'\1_\2\t\3\t\4'
               cmd = re.compile(search_str)
               sys.stdout.write('#name\tlength\tcoverage\n')
               for i,line in enumerate(sys.stdin):
                  if cmd.match(line):
                     sys.stdout.write(cmd.sub(replace_str,line))

####################
      - entryname: spades_wrapper.sh
        entry: |
               ####################spades launcher
               #!/bin/bash
               ## An example command looks like:
               ## spades.py -k 21,33,55,77,99,127 --careful -1 Y.fastq.gz -2 X.fastq.gz -t 24 -o output 
               read -d '' MEMORY_GB << EOF
               ${ 
                //compute memory limits
                var mem="250"
                if(runtime.ram){
                  var bt=runtime.ram * 1048576 //mebibytes to bytes               
                  bt=bt*0.000000001 //bytes to Gigabytes 
                  mem=parseInt(bt)
                }
                return "-"+mem
               }
               EOF
               #echo "MEMORY_GB limit: \$MEMORY_GB"
               read -d '' CORES << EOF
               ${ 
                //compute core / slot
                var cores="16" 
                if(runtime.cores){
                  cores=runtime.cores
                }
                return  cores
               }
               EOF
               #echo "CORES : \$CORES"
               read -d '' SPADES_OPT << EOF
               ${   
                 var opt=""
                   if (inputs.sc==true){
                    opt+=" --sc "
                   }
                   if (inputs.onlyassembler==true){
                    opt+=" --only-assembler "
                   }
                   if (inputs.careful==true){
                    opt+=" --careful "
                   }
                 return opt  
               }
               EOF


               read -d '' CMD_BASE << EOF
               spades.py -o . --disable-gzip-output  
                 \$SPADES_OPT 
                -t \$CORES -m \$MEMORY_GB
               ${   
                 var opt=""
                   if (inputs.auto_kmer_choice==false){
                    opt+=" -k "+inputs.kmers 
                   }
                 return opt  
               }
               ${
                 var opt=""
                 if (inputs.cov_state==null || inputs.cov_state == "auto"){
                      opt+=" --cov-cutoff 'auto' "
                  }
                 else if (inputs.cov_state == "value"){
                      opt+=" --cov-cutoff '"+inputs.cov_cutoff+"' "
                  }
                 if (inputs.iontorrent==true ){
                      opt+=" --iontorrent "
                 } 
                 return opt  
               }
               EOF
               ##########################################
               ##Sequence files from libraries
               read -d '' CMD_READ1 << EOF
               ${
               var opt=""
               var lib_prefix = {}
               var lib_meta=inputs.libraries_metadata   
               if (lib_meta!=null){
                   for(var j=0; j<lib_meta.length;j++){
                      var prefix=""
                      var lmeta=lib_meta[j]
                      if (lmeta.lib_type !=null && lmeta.lib_type == "paired_end"){
                           prefix = 'pe'
                       } 
                      else if (lmeta.lib_type !=null && lmeta.lib_type == "mate_paired"){
                           prefix = 'mp'
                       }
                      else if (lmeta.lib_type !=null && lmeta.lib_type == "nxmate_paired"){
                           prefix = 'nxmate'
                       }                       
                      else {
                           prefix = 'hqmp'
                       }
                      var idx=lmeta.lib_index
                      lib_prefix[idx]=prefix
                      opt+=" --"+prefix+idx+"-"+lmeta.orientation+" "
                  } 
                }  
               var libraries = []
               if(inputs.libraries_fwd_rev!=null){ 
                  for (var i = 0; i < inputs.libraries_fwd_rev.length ; i++) {
                     var lib=inputs.libraries_fwd_rev[i] 
                     lib.file_type="separate"
                     libraries[i] = lib
                  }
                } 
               if(inputs.libraries_mono!=null){
                  for (var i = 0; i < inputs.libraries_mono.length ; i++) {
                    var ei= libraries.length
                    libraries[ei] = inputs.libraries_mono[i]
                  }
               }  
               for(var j=0; j<libraries.length;j++){
                 var lib=libraries[j]
                 var idx=lib.lib_index  
                 var prefix=lib_prefix[idx]
                 if(lib.file_type!=null){
                    if ( lib.file_type == "separate"){
                      opt+=" --"+prefix+idx+"-1 "+"fastq:"+lib.fwd_reads.path+" "
                      opt+=" --"+prefix+idx+"-2 "+"fastq:"+lib.rev_reads.path+" "
                    }else{
                      var suffix
                      if ( lib.file_type == "interleaved"){
                        suffix="12"
                      }
                      if ( lib.file_type == "merged"){
                        suffix="m"
                      }
                      if ( lib.file_type == "unpaired"){
                        suffix="s"
                      }
                      opt+=" --"+prefix+idx+"-"+suffix+" "+"fastq:"+lib.reads.path+" "
                    }
                  }
                }  
                 return opt  
               }
               EOF
               ##########################################



               #########################>DEBUG

               ##########################################
               ##Sequence files from libraries
               read -d '' ZZDEBUG << EOF
               ${
               //var opt="inputs.libraries_metadata.length:"+inputs.libraries_metadata.length+" "              
               //var str = JSON.stringify(inputs.libraries_metadata)
               var opt=""
               var lib_prefix = {}
               var lib_meta=inputs.libraries_metadata   
               if (lib_meta!=null){
                   opt+=" not null "
                   for(var j=0; j<lib_meta.length;j++){
                      var prefix=""
                      var lmeta=lib_meta[j]
 
                      opt+=" !!-AA-"+lmeta.lib_index+"-"+lmeta.orientation+" "
                  } 
                }  
               var libraries = []
               if(inputs.libraries_fwd_rev!=null){ 
                  for (var i = 0; i < inputs.libraries_fwd_rev.length ; i++) {
                     var lib=inputs.libraries_fwd_rev[i] 
                     lib.file_type="separate"
                     libraries[i] = lib
                  }
                } 
               if(inputs.libraries_mono!=null){
                  for (var i = 0; i < inputs.libraries_mono.length ; i++) {
                    var ei= libraries.length
                    libraries[ei] = inputs.libraries_mono[i]
                  }
               }  
               for(var j=0; j<libraries.length;j++){
                 var lib=libraries[j]
                 var idx=lib.lib_index  
                 var prefix=lib_prefix[idx]
                 if(lib.file_type!=null){
                    if ( lib.file_type == "separate"){
                      opt+=" --"+prefix+idx+"-1 "+"fastq:"+lib.fwd_reads.path+" "
                      opt+=" --"+prefix+idx+"-2 "+"fastq:"+lib.rev_reads.path+" "
                    }else{
                      var suffix
                      if ( lib.file_type == "interleaved"){
                        suffix="12"
                      }
                      if ( lib.file_type == "merged"){
                        suffix="m"
                      }
                      if ( lib.file_type == "unpaired"){
                        suffix="s"
                      }
                      opt+=" --"+prefix+idx+"-"+suffix+" "+"fastq:"+lib.reads.path+" "
                    }
                  }
                }  
                 return opt  
               }
               EOF
               ##########################################


               #########################<DEBUG










               read -d '' CMD_READ2 << EOF
               ${
                 var opt=""
                 if (inputs.pacbio_reads!=null){
                   for(var i=0; i<inputs.pacbio_reads.length;i++){
                      var read=inputs.pacbio_reads[i]
                      opt+=" --pacbio fastq:"+read.path+" "
                   }
                  }
                 if (inputs.nanopore_reads!=null){
                   for(var i=0; i<inputs.nanopore_reads.length;i++){
                      var read=inputs.nanopore_reads[i]
                      opt+=" --nanopore fastq:"+read.path+" "
                   }
                  }
                 if (inputs.sanger_reads!=null){
                   for(var i=0; i<inputs.sanger_reads.length;i++){
                      var read=inputs.sanger_reads[i]
                      opt+=" --sanger fastq:"+read.path+" "
                   }
                  }
                 if (inputs.trusted_contigs!=null){
                   for(var i=0; i<inputs.trusted_contigs.length;i++){
                      var read=inputs.trusted_contigs[i]
                      opt+=" --trusted-contigs fastq:"+read.path+" "
                   }
                  }
                 if (inputs.untrusted_contigs!=null){
                   for(var i=0; i<inputs.untrusted_contigs.length;i++){
                      var read=inputs.untrusted_contigs[i]
                      opt+=" --untrusted-contigs fastq:"+read.path+" "
                   }
                  }
                 return opt  
               }
               EOF
               ##########################################
               read -d '' CMD_POST << EOF
                 && python write_tsv.py < contigs.fasta > out_contig_stats.tab
                 && python write_tsv.py < scaffolds.fasta > out_scaffold_stats.tab
               EOF
               ##########################################   
               # echo "CMD_BASE : \$CMD_BASE"  > zz.txt   \\
               #        && echo "CMD_READ1 : \$CMD_READ1"  >> zz.txt   \\
               #        && echo "CMD_READ2 : \$CMD_READ2"  >> zz.txt   \\
               #        && echo "CMD_POST : \$CMD_POST"   >> zz.txt
               COMMAND="\$CMD_BASE \$CMD_READ1 \$CMD_READ2 \$CMD_POST"
               COMMAND=\$(echo $COMMAND|tr -d '\\n')
               echo "\$COMMAND"  > run_spades.sh
               bash  ./run_spades.sh

 
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biocontainers/unicycler:v0.4.7dfsg-2-deb_cv1'
  
doc: |
  example workflow for js wrapper generation
  see  https://github.com/rrwick/Unicycler
  outputs : all genretated files
