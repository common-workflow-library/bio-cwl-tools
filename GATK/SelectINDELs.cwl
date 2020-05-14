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
            "id": "mutantVCF"
        }, 
        {
            "type": "File", 
            "id": "ReferenceGenome"
        }, 
        {
            "type": "File", 
            "id": "RefIndex"
        }, 
        {
            "type": "File", 
            "id": "RefDict"
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
                "glob": "$(inputs.sampleName).mutantindel.vcf"
            }, 
            "type": "File", 
            "id": "rawVCF"
        }
    ], 
    "baseCommand": [], 
    "class": "CommandLineTool", 
    "arguments": [
        {
            "shellQuote": false, 
            "valueFrom": "$(inputs.GATK.path) SelectVariants -R $(inputs.ReferenceGenome.path) -V $(inputs.mutantVCF.path) -O $(inputs.sampleName).mutantindel.vcf -select-type-to-include INDEL"
        }
    ], 
    "id": "SelectINDELs"
}