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
                "glob": "$(inputs.sampleName).mutantfilter.vcf"
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
            "valueFrom": "$(inputs.GATK.path) IndexFeatureFile -F $(inputs.mutantVCF.path)\r\r $(inputs.GATK.path) VariantFiltration -R $(inputs.ReferenceGenome.path) -V $(inputs.mutantVCF.path) -window 35 -cluster 3 -O $(inputs.sampleName).mutantfilter.vcf"
        }
    ], 
    "id": "VariantFilteration"
}