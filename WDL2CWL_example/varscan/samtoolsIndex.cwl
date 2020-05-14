#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).bam.bai"
            },
            "type": "File",
            "id": "rawBAM"
        }
    ],
    "inputs": [
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
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "samtools index $(inputs.inputBAM.path) > $(inputs.sampleName).bam.bai"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "samtoolsIndex"
}