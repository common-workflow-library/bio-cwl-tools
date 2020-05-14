#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).sorted.bam"
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
        },
        {
            "type": "string",
            "id": "docker"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "dockerPull": "docker",
            "class": "DockerRequirement"
        }
    ],
    "arguments": [
        {
            "shellQuote": false,
            "valueFrom": "samtools sort -l 0 -o $(inputs.sampleName).sorted.bam $(inputs.inputBAM.path)"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "samtoolsSort"
}