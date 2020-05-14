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
        }
    ],
    "class": "CommandLineTool",
    "inputs": [
        {
            "id": "inputSAM",
            "type": "File"
        },
        {
            "id": "sampleName",
            "type": "string"
        },
        {
            "id": "samtoolsPath",
            "type": "string"
        }
    ],
    "id": "samtoolsView",
    "outputs": [
        {
            "id": "rawBAM",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.sampleName).bam"
            }
        }
    ],
    "arguments": [
        {
            "valueFrom": "samtools view -bS $(inputs.inputSAM.path) > $(inputs.sampleName).bam",
            "shellQuote": false
        }
    ],
    "cwlVersion": "v1.0"
}