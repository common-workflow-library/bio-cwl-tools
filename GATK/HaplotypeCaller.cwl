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
                "glob": "$(inputs.sampleName).mutant.vcf"
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
            "valueFrom": "samtools index $(inputs.inputBAM.path) > $(inputs.sampleName).split.bam.bai\r\r                      \r\r                $(inputs.GATK.path) HaplotypeCaller -R $(inputs.ReferenceGenome.path) -I $(inputs.inputBAM.path) -O $(inputs.sampleName).mutant.vcf"
        }
    ], 
    "id": "HaplotypeCaller"
}