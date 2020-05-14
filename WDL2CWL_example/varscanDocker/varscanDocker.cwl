#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "steps": [
        {
            "in": [
                {
                    "source": "referenceGenome",
                    "id": "ReferenceGenome"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                },
                {
                    "source": "bowtieDocker",
                    "id": "docker"
                },
                {
                    "source": "name",
                    "id": "index"
                },
                {
                    "source": "Alignment_leftFastq",
                    "id": "leftFastq"
                },
                {
                    "source": "Alignment_rightFastq",
                    "id": "rightFastq"
                }
            ],
            "run": "Alignment.cwl",
            "out": [
                {
                    "id": "rawSAM"
                }
            ],
            "id": "Alignment"
        },
        {
            "in": [
                {
                    "source": "#Alignment/rawSAM",
                    "id": "inputSAM"
                },
                {
                    "source": "samtoolsDocker",
                    "id": "docker"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                }
            ],
            "run": "samtoolsView.cwl",
            "out": [
                {
                    "id": "rawBAM"
                }
            ],
            "id": "samtoolsView"
        },
        {
            "in": [
                {
                    "source": "#samtoolsView/rawBAM",
                    "id": "inputBAM"
                },
                {
                    "source": "samtoolsDocker",
                    "id": "docker"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                }
            ],
            "run": "samtoolsSort.cwl",
            "out": [
                {
                    "id": "rawBAM"
                }
            ],
            "id": "samtoolsSort"
        },
        {
            "in": [
                {
                    "source": "#samtoolsSort/rawBAM",
                    "id": "inputBAM"
                },
                {
                    "source": "samtoolsDocker",
                    "id": "docker"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                }
            ],
            "run": "samtoolsIndex.cwl",
            "out": [
                {
                    "id": "rawBAM"
                }
            ],
            "id": "samtoolsIndex"
        },
        {
            "in": [
                {
                    "source": "referenceGenome",
                    "id": "ReferenceGenome"
                },
                {
                    "source": "#samtoolsSort/rawBAM",
                    "id": "inputBAM"
                },
                {
                    "source": "samtoolsDocker",
                    "id": "docker"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                }
            ],
            "run": "samtoolsMpileup.cwl",
            "out": [
                {
                    "id": "rawMpileup"
                }
            ],
            "id": "samtoolsMpileup"
        },
        {
            "in": [
                {
                    "source": "#samtoolsMpileup/rawMpileup",
                    "id": "inputMpileup"
                },
                {
                    "source": "name",
                    "id": "sampleName"
                },
                {
                    "source": "varscanDocker",
                    "id": "docker"
                }
            ],
            "run": "mpileup2snp.cwl",
            "out": [
                {
                    "id": "snpVCF"
                }
            ],
            "id": "mpileup2snp"
        }
    ],
    "outputs": [
        {
            "outputSource": "#Alignment/rawSAM",
            "type": "File",
            "id": "Alignment_rawSAM"
        },
        {
            "outputSource": "#samtoolsView/rawBAM",
            "type": "File",
            "id": "samtoolsView_rawBAM"
        },
        {
            "outputSource": "#samtoolsSort/rawBAM",
            "type": "File",
            "id": "samtoolsSort_rawBAM"
        },
        {
            "outputSource": "#samtoolsIndex/rawBAM",
            "type": "File",
            "id": "samtoolsIndex_rawBAM"
        },
        {
            "outputSource": "#samtoolsMpileup/rawMpileup",
            "type": "File",
            "id": "samtoolsMpileup_rawMpileup"
        },
        {
            "outputSource": "#mpileup2snp/snpVCF",
            "type": "File",
            "id": "mpileup2snp_snpVCF"
        }
    ],
    "inputs": [
        {
            "type": "string",
            "id": "varscanDocker"
        },
        {
            "type": "string",
            "id": "bowtieDocker"
        },
        {
            "type": "string",
            "id": "samtoolsDocker"
        },
        {
            "type": "File",
            "id": "referenceGenome"
        },
        {
            "type": "string",
            "id": "name"
        },
        {
            "type": "File",
            "id": "Alignment_leftFastq"
        },
        {
            "type": "File",
            "id": "Alignment_rightFastq"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "cwlVersion": "v1.0",
    "class": "Workflow",
    "id": "varscanDocker"
}