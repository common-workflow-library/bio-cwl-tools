#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).mpileup"
            },
            "type": "File",
            "id": "rawMpileup"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "inputBAM"
        },
        {
            "type": "File",
            "id": "ReferenceGenome"
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
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "samtools index $(inputs.inputBAM.path)samtools mpileup -B -f $(inputs.ReferenceGenome.path) $(inputs.inputBAM.path) > $(inputs.sampleName).mpileup"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "samtoolsMpileup"
}