#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "docker"
        }
    ],
    "class": "CommandLineTool",
    "inputs": [
        {
            "id": "inputBAM",
            "type": "File"
        },
        {
            "id": "sampleName",
            "type": "string"
        },
        {
            "id": "ReferenceGenome",
            "type": "File"
        },
        {
            "id": "docker",
            "type": "string"
        }
    ],
    "id": "bcftools",
    "outputs": [
        {
            "id": "rawVCF",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.sampleName).vcf"
            }
        }
    ],
    "arguments": [
        {
            "valueFrom": "bcftools mpileup -f $(inputs.ReferenceGenome.path) $(inputs.inputBAM.path) | bcftools call -mv -Ob -o $(inputs.sampleName).vcf",
            "shellQuote": false
        }
    ],
    "cwlVersion": "v1.0"
}