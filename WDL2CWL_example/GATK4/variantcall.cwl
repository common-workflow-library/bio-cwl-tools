#!/usr/bin/env cwl-runner
# This tool description was generated automatically by wdl2cwl ver. 0.2

{
    "cwlVersion": "v1.0", 
    "inputs": [
        {
            "type": "File", 
            "id": "gatk"
        }, 
        {
            "type": "File", 
            "id": "refIndex"
        }, 
        {
            "type": "File", 
            "id": "refDict"
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
    "outputs": [
        {
            "outputSource": "#Alignment/rawSAM", 
            "type": "File", 
            "id": "Alignment_rawSAM"
        }, 
        {
            "outputSource": "#AddOrReplaceReadGroups/rawBAM", 
            "type": "File", 
            "id": "AddOrReplaceReadGroups_rawBAM"
        }, 
        {
            "outputSource": "#SortSam/rawBAM", 
            "type": "File", 
            "id": "SortSam_rawBAM"
        }, 
        {
            "outputSource": "#ReferenceSeqIndex/refIndex", 
            "type": "File", 
            "id": "ReferenceSeqIndex_refIndex"
        }, 
        {
            "outputSource": "#ReferenceSeqDictionary/refDict", 
            "type": "File", 
            "id": "ReferenceSeqDictionary_refDict"
        }, 
        {
            "outputSource": "#MarkDuplicates/rawBAM", 
            "type": "File", 
            "id": "MarkDuplicates_rawBAM"
        }, 
        {
            "outputSource": "#SplitNCigarReads/rawBAM", 
            "type": "File", 
            "id": "SplitNCigarReads_rawBAM"
        }, 
        {
            "outputSource": "#HaplotypeCaller/rawVCF", 
            "type": "File", 
            "id": "HaplotypeCaller_rawVCF"
        }, 
        {
            "outputSource": "#VariantFilteration/rawVCF", 
            "type": "File", 
            "id": "VariantFilteration_rawVCF"
        }, 
        {
            "outputSource": "#SelectSNPs/rawVCF", 
            "type": "File", 
            "id": "SelectSNPs_rawVCF"
        }, 
        {
            "outputSource": "#SelectINDELs/rawVCF", 
            "type": "File", 
            "id": "SelectINDELs_rawVCF"
        }
    ], 
    "class": "Workflow", 
    "steps": [
        {
            "out": [
                {
                    "id": "rawSAM"
                }
            ], 
            "run": "Alignment.cwl", 
            "id": "Alignment", 
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
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawBAM"
                }
            ], 
            "run": "AddOrReplaceReadGroups.cwl", 
            "id": "AddOrReplaceReadGroups", 
            "in": [
                {
                    "source": "#Alignment/rawSAM", 
                    "id": "inputSAM"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawBAM"
                }
            ], 
            "run": "SortSam.cwl", 
            "id": "SortSam", 
            "in": [
                {
                    "source": "#AddOrReplaceReadGroups/rawBAM", 
                    "id": "inputBAM"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "refIndex"
                }
            ], 
            "run": "ReferenceSeqIndex.cwl", 
            "id": "ReferenceSeqIndex", 
            "in": [
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "refDict"
                }
            ], 
            "run": "ReferenceSeqDictionary.cwl", 
            "id": "ReferenceSeqDictionary", 
            "in": [
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawBAM"
                }
            ], 
            "run": "MarkDuplicates.cwl", 
            "id": "MarkDuplicates", 
            "in": [
                {
                    "source": "#SortSam/rawBAM", 
                    "id": "inputBAM"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawBAM"
                }
            ], 
            "run": "SplitNCigarReads.cwl", 
            "id": "SplitNCigarReads", 
            "in": [
                {
                    "source": "#MarkDuplicates/rawBAM", 
                    "id": "inputBAM"
                }, 
                {
                    "source": "refIndex", 
                    "id": "RefIndex"
                }, 
                {
                    "source": "refDict", 
                    "id": "RefDict"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawVCF"
                }
            ], 
            "run": "HaplotypeCaller.cwl", 
            "id": "HaplotypeCaller", 
            "in": [
                {
                    "source": "#SplitNCigarReads/rawBAM", 
                    "id": "inputBAM"
                }, 
                {
                    "source": "refIndex", 
                    "id": "RefIndex"
                }, 
                {
                    "source": "refDict", 
                    "id": "RefDict"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawVCF"
                }
            ], 
            "run": "VariantFilteration.cwl", 
            "id": "VariantFilteration", 
            "in": [
                {
                    "source": "#HaplotypeCaller/rawVCF", 
                    "id": "mutantVCF"
                }, 
                {
                    "source": "refIndex", 
                    "id": "RefIndex"
                }, 
                {
                    "source": "refDict", 
                    "id": "RefDict"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawVCF"
                }
            ], 
            "run": "SelectSNPs.cwl", 
            "id": "SelectSNPs", 
            "in": [
                {
                    "source": "#VariantFilteration/rawVCF", 
                    "id": "mutantVCF"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "refIndex", 
                    "id": "RefIndex"
                }, 
                {
                    "source": "refDict", 
                    "id": "RefDict"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }, 
        {
            "out": [
                {
                    "id": "rawVCF"
                }
            ], 
            "run": "SelectINDELs.cwl", 
            "id": "SelectINDELs", 
            "in": [
                {
                    "source": "#VariantFilteration/rawVCF", 
                    "id": "mutantVCF"
                }, 
                {
                    "source": "gatk", 
                    "id": "GATK"
                }, 
                {
                    "source": "referenceGenome", 
                    "id": "ReferenceGenome"
                }, 
                {
                    "source": "refIndex", 
                    "id": "RefIndex"
                }, 
                {
                    "source": "refDict", 
                    "id": "RefDict"
                }, 
                {
                    "source": "name", 
                    "id": "sampleName"
                }
            ]
        }
    ], 
    "id": "variantcall"
}