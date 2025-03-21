#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Kallisto quant
doc: |

  Docs: https://pachterlab.github.io/kallisto/

  Computes equivalence classes for reads and quantifies abundances

  Usage: kallisto quant [arguments] FASTQ-files

  Required arguments:
  -i, --index=STRING            Filename for the kallisto index to be used for
                                quantification
  -o, --output-dir=STRING       Directory to write output to

  Optional arguments:
  -b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
      --seed=INT                Seed for the bootstrap sampling (default: 42)
      --plaintext               Output plaintext instead of HDF5
      --single                  Quantify single-end reads
      --single-overhang         Include reads where unobserved rest of fragment is
                                predicted to lie outside a transcript
      --fr-stranded             Strand specific reads, first read forward
      --rf-stranded             Strand specific reads, first read reverse
  -l, --fragment-length=DOUBLE  Estimated average fragment length
  -s, --sd=DOUBLE               Estimated standard deviation of fragment length
                                (default: -l, -s values are estimated from paired
                                end data, but are required when using --single)
  -p, --priors                  Priors for the EM algorithm, either as raw counts or as
                                probabilities. Pseudocounts are added to raw reads to
                                prevent zero valued priors. Supplied in the same order
                                as the transcripts in the transcriptome
  -t, --threads=INT             Number of threads to use (default: 1)
      --verbose                 Print out progress information every 1M proccessed reads


  This CWL was adapted from: https://github.com/common-workflow-library/bio-cwl-tools/commit/91c42fb809ce18eafe16155cca0abf362270c0fe


hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1
  SoftwareRequirement:
    packages:
      - package: kallisto
        version: [ "0.51.1" ]
        specs: [ https://identifiers.org/biotools/kallisto ]

inputs:
  InputReads:
    type: File[]
    format: edam:format_1930  # FASTQ
    inputBinding:
      position: 200

  QuantOutfolder: 
    type: string

  Index:
    type: File
    inputBinding:
      position: 1
      prefix: "--index"

  isSingle:
    type: boolean
    inputBinding:
      position: 2
      prefix: "--single"

  #Optional Inputs

  isBias:
    type: boolean?
    inputBinding:
      prefix: "--bias"

  isFusion:
    type: boolean?
    inputBinding:
      prefix: "--fusion"

  isSingleOverhang:
    type: boolean?
    inputBinding:
      prefix: "--single-overhang"
  
  FragmentLength:
    type: double?
    inputBinding:
      separate: false
      prefix: "--fragment-length="
  
  StandardDeviation:
    type: double?
    inputBinding:
      prefix: "--sd"
  
  BootstrapSamples:
    type: int?
    inputBinding:
      separate: false
      prefix: "--bootstrap-samples="
  
  Seed:
    type: int?
    inputBinding:
      prefix: "--seed"

#Using record inputs to create mutually exclusive inputs
  Strand:
    type:
      - "null"
      - type: record
        name: forward
        fields:
          forward:
              type: boolean
              inputBinding:
                prefix: "--fr-stranded"

      - type: record
        name: reverse
        fields:
          reverse:
            type: boolean
            inputBinding:
              prefix: "--rf-stranded"

  PseudoBam:
    type: boolean?
    inputBinding:
      prefix: "--pseudobam"

#Using record inputs to create dependent inputs
  
  GenomeBam:
    type:
      - "null"
      - type: record
        name: genome_bam
        fields:
          genomebam:
            type: boolean
            inputBinding:
              prefix: "--genomebam"

          gtf:
            type: File
            inputBinding:
              prefix: "--gtf"

          chromosomes:
            type: File
            inputBinding:
              prefix: "--chromosomes"

baseCommand: [ kallisto, quant ]

arguments: [ "--output-dir", $(inputs.QuantOutfolder) ]

outputs:

  kallistoQuantOutDir:
    type: Directory
    outputBinding:
      glob: $(runtime.outdir)/$(inputs.QuantOutfolder)


$namespaces:
  edam: https://edamontology.org/
  s: https://schema.org/
$schemas:
  - https://edamontology.org/EDAM_1.25.owl
  - https://schema.org/version/latest/schemaorg-current-https.rdf

s:license: https://spdx.org/licenses/BSD-2-Clause
s:citation: https://dx.doi.org/10.1038/nbt.3519
s:codeRepository: https://github.com/pachterlab/kallisto
