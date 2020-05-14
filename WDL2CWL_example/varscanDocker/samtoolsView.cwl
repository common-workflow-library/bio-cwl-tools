#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).bam"
            },
            "type": "File",
            "id": "rawBAM"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "inputSAM"
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
            "valueFrom": "samtools view -bS $(inputs.inputSAM.path) > $(inputs.sampleName).bam"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "samtoolsView"
}