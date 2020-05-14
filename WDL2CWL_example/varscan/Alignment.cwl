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
            "valueFrom": "export PATH=$PATH:/home/ngsap2/Downloads/bowtie2-2.4.1          bowtie2-build $(inputs.ReferenceGenome.path) $(inputs.index)                     bowtie2 -q -x $(inputs.index) -1 $(inputs.leftFastq.path) -2 $(inputs.rightFastq.path) -S $(inputs.sampleName).sam"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "CommandLineTool",
    "id": "Alignment"
}