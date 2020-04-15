#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: "viterbi: Viterbi realignment
Probabilistic realignment of your already mapped reads, which corrects mapping errors (run right after mapping). Not recommended for non-Illumina data."

hints:
  RessourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: quay.io/biocontainers/lofreq:2.1.4--py27hc3dfafe_1

requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference)
    
baseCommand: [lofreq, viterbi]

arguments:
  - prefix: --out
    valueFrom: $(inputs.reads.nameroot)_realigned.bam
    position: 99

inputs:
  reference:
    type: File
  #  format:
  #    - edam:format_ #FASTA
    inputBinding:
      prefix: --ref
      position: 1
      
  reads:  
    type: File
  #  format:
  #    - edam:format_2572
    inputBinding:
      position: 100

  keepflags:
    type: boolean?
    doc: "Delete flags MC, MD, NM, and A? These flags are all prone to getting invalidated during realignment. Keep them only if you know what you are doing."
    inputBinding:
      prefix: --keepflags
      position: 2   
    default: false
      
  bq2_handling:
    type:
      - 'null'
      - type: enum
        symbols:
          - keep
          - dynamic
          - fixed
    inputBinding:
      position: 3
    doc: "How to handle base qualities of 2? In sequenced reads obtained with Illumina sequencing pipelines before version 1.8, base quality 2 is special in that it serves as a general indicator of low quality of the corresponding bases. For such reads, the tool can make an optimistic guess of the real quality of such bases by replacing base qualities of 2 with the median of all other base qualities observed in the read. Alternatively, you can provide a fixed replacement value. For recently obtained sequencing data, just keep BQ2 values unchanged (the default) since they have no special meaning."

  defqual:
    doc: "If bq2_handling=fixed"
    type: int?
    inputBinding:
      prefix: --defqual
      position: 4  
outputs:
  realigned:
    type: File
  #  format:
  #    - edam:format_2572

    outputBinding:
      glob: $(inputs.reads.nameroot)_realigned.bam


#  ${adv_options.delflags}
#    --defqual ${adv_options.bq2_handling.defqual}


#    samtools sort -T "\${TMPDIR:-.}" -@ \${GALAXY_SLOTS:-1} -O BAM -o '$realigned' tmp.bam

#$namespaces:
#  edam: http://edamontology.org/
#$schemas:
#  - http://edamontology.org/EDAM_1.18.owl

