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
            "id": "referenceGenome",
            "type": "File"
        },
        {
            "id": "name",
            "type": "string"
        },
        {
            "id": "samtoolsPath",
            "type": "string"
        },
        {
            "id": "bcftoolsPath",
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
    "id": "samtools",
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
                    "id": "samtoolsPath",
                    "source": "samtoolsPath"
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
                    "id": "samtoolsPath",
                    "source": "samtoolsPath"
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
                    "id": "samtoolsPath",
                    "source": "samtoolsPath"
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
                    "id": "bcftoolsPath",
                    "source": "bcftoolsPath"
                },
                {
                    "id": "ReferenceGenome",
                    "source": "referenceGenome"
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