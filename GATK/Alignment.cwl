#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "cwlVersion": "v1.0", 
    "inputs": [
        {
            "type": "File", 
            "id": "leftFastq"
        }, 
        {
            "type": "File", 
            "id": "rightFastq"
        }, 
        {
            "type": "File", 
            "id": "ReferenceGenome"
        }, 
        {
            "type": "string", 
            "id": "sampleName"
        }, 
        {
            "type": "string", 
            "id": "index"
        }
    ], 
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        }, 
        {
            "class": "InlineJavascriptRequirement"
        }
    ], 
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).sam"
            }, 
            "type": "File", 
            "id": "rawSAM"
        }
    ], 
    "baseCommand": [], 
    "class": "CommandLineTool", 
    "arguments": [
        {
            "shellQuote": false, 
            "valueFrom": "export PATH=$PATH:/home/../../bowtie\r\r          bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index) \r\r          \r\r          bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam"
        }
    ], 
    "id": "Alignment"
}