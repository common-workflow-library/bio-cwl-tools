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
            "id": "inputBAM",
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
    "id": "samtoolsSort",
    "outputs": [
        {
            "id": "rawBAM",
            "type": "File",
            "outputBinding": {
                "glob": "$(inputs.sampleName).sorted.bam"
            }
        }
    ],
    "arguments": [
        {
            "valueFrom": "samtools sort -l 0 -o $(inputs.sampleName).sorted.bam $(inputs.inputBAM.path)",
            "shellQuote": false
        }
    ],
    "cwlVersion": "v1.0"
}