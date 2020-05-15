#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

    cwlVersion: v1.0
    baseCommand: []
    requirements:         
           - class: ShellCommandRequirement
        
        
           - class: InlineJavascriptRequirement
        
        
           - class: DockerRequirement
             dockerPull: docker
        
    
    class: CommandLineTool
    inputs: 
        leftFastq:
            id: leftFastq
            type: File
        
        rightFastq:
            id: rightFastq
            type: File
        
        ReferenceGenome:
             id: ReferenceGenome
             type: File
        
        sampleName:
             id: sampleName
             type: string
        
        index:
             id: index
             type: string
        
        docker:
             id: docker
             type: string
        
    
    id: Alignment
    outputs:
        
             id: rawSAM
             type: File
             outputBinding:
                glob: $(inputs.sampleName).sam
      
    arguments: 
             valueFrom: bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index) | bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam
             shellQuote: false
        
    
    
