#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "outputs": [
        {
            "id": "Alignment_rawSAM",
            "type": "File",
            "outputSource": "#Alignment/rawSAM"
        },
        {
            "id": "samtoolsView_rawBAM",
            "type": "File",
            "outputSource": "#samtoolsView/rawBAM"
        },
        {
            "id": "samtoolsSort_rawBAM",
            "type": "File",
            "outputSource": "#samtoolsSort/rawBAM"
        },
        {
            "id": "samtoolsIndex_rawBAM",
            "type": "File",
            "outputSource": "#samtoolsIndex/rawBAM"
        },
        {
            "id": "bcftools_rawVCF",
            "type": "File",
            "outputSource": "#bcftools/rawVCF"
        }
    ],
    "class": "Workflow",
    "inputs": [
        {
            "id": "bowtieDocker",
            "type": "string"
        },
        {
            "id": "samtoolsDocker",
            "type": "string"
        },
        {
            "id": "bcftoolsDocker",
            "type": "string"
        },
        {
            "id": "referenceGenome",
            "type": "File"
        },
        {
            "id": "name",
            "type": "string"
        },
        {
            "id": "Alignment_leftFastq",
            "type": "File"
        },
        {
            "id": "Alignment_rightFastq",
            "type": "File"
        }
    ],
    "id": "samtoolsDocker",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "steps": [
        {
            "id": "Alignment",
            "run": "Alignment.cwl",
            "in": [
                {
                    "id": "ReferenceGenome",
                    "source": "referenceGenome"
                },
                {
                    "id": "sampleName",
                    "source": "name"
                },
                {
                    "id": "docker",
                    "source": "bowtieDocker"
                },
                {
                    "id": "index",
                    "source": "name"
                },
                {
                    "id": "leftFastq",
                    "source": "Alignment_leftFastq"
                },
                {
                    "id": "rightFastq",
                    "source": "Alignment_rightFastq"
                }
            ],
            "out": [
                {
                    "id": "rawSAM"
                }
            ]
        },
        {
            "id": "samtoolsView",
            "run": "samtoolsView.cwl",
            "in": [
                {
                    "id": "inputSAM",
                    "source": "#Alignment/rawSAM"
                },
                {
                    "id": "docker",
                    "source": "samtoolsDocker"
                },
                {
                    "id": "sampleName",
                    "source": "name"
                }
            ],
            "out": [
                {
                    "id": "rawBAM"
                }
            ]
        },
        {
            "id": "samtoolsSort",
            "run": "samtoolsSort.cwl",
            "in": [
                {
                    "id": "inputBAM",
                    "source": "#samtoolsView/rawBAM"
                },
                {
                    "id": "docker",
                    "source": "samtoolsDocker"
                },
                {
                    "id": "sampleName",
                    "source": "name"
                }
            ],
            "out": [
                {
                    "id": "rawBAM"
                }
            ]
        },
        {
            "id": "samtoolsIndex",
            "run": "samtoolsIndex.cwl",
            "in": [
                {
                    "id": "inputBAM",
                    "source": "#samtoolsSort/rawBAM"
                },
                {
                    "id": "docker",
                    "source": "samtoolsDocker"
                },
                {
                    "id": "sampleName",
                    "source": "name"
                }
            ],
            "out": [
                {
                    "id": "rawBAM"
                }
            ]
        },
        {
            "id": "bcftools",
            "run": "bcftools.cwl",
            "in": [
                {
                    "id": "inputBAM",
                    "source": "#samtoolsSort/rawBAM"
                },
                {
                    "id": "sampleName",
                    "source": "name"
                },
                {
                    "id": "ReferenceGenome",
                    "source": "referenceGenome"
                },
                {
                    "id": "docker",
                    "source": "bcftoolsDocker"
                }
            ],
            "out": [
                {
                    "id": "rawVCF"
                }
            ]
        }
    ],
    "cwlVersion": "v1.0"
}