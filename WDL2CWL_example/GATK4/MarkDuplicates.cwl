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
            "id": "inputBAM"
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
                "glob": "$(inputs.sampleName).markdup.bam"
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
            "valueFrom": "$(inputs.GATK.path) MarkDuplicates -I $(inputs.inputBAM.path) -O $(inputs.sampleName).markdup.bam -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT -M output.metrics"
        }
    ], 
    "id": "MarkDuplicates"
}