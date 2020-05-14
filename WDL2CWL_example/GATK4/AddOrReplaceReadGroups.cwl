#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "cwlVersion": "v1.0", 
    "inputs": [
        {
            "type": "File", 
            "id": "GATK"
        }, 
        {
            "type": "File", 
            "id": "inputSAM"
        }, 
        {
            "type": "string", 
            "id": "sampleName"
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
                "glob": "$(inputs.sampleName).bam"
            }, 
            "type": "File", 
            "id": "rawBAM"
        }
    ], 
    "baseCommand": [], 
    "class": "CommandLineTool", 
    "arguments": [
        {
            "shellQuote": false, 
            "valueFrom": "$(inputs.GATK.path) AddOrReplaceReadGroups -I $(inputs.inputSAM.path) -O $(inputs.sampleName).bam -RGID 1 -RGLB 445_LIB -RGPL illumina -RGSM RNA -RGPU illumina"
        }
    ], 
    "id": "AddOrReplaceReadGroups"
}