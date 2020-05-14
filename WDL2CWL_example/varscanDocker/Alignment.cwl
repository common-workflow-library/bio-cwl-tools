#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "baseCommand": [],
    "outputs": [
        {
            "outputBinding": {
                "glob": "$(inputs.sampleName).sam"
            },
            "type": "File",
            "id": "rawSAM"
        }
    ],
    "inputs": [
        {
            "type": "File",
            "id": "leftFastq"
        },
        {
            "type": "File",
            "id": "rightFastq"
        },
        {
            "type": "File",
            "id": "ReferenceGenome"
        },
        {
            "type": "string",
            "id": "sampleName"
        },
        {
            "type": "string",
            "id": "index"
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
            "valueFrom": "bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index)                     bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "Alignment"
}