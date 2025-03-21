#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: Kallisto quant
doc: |2

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

inputs:
  BootstrapSamples:
    type: int?
    inputBinding:
      prefix: --bootstrap-samples=
      separate: false
  FragmentLength:
    type: double?
    inputBinding:
      prefix: --fragment-length=
      separate: false
  GenomeBam:
    type:
    - 'null'
    - name: genome_bam
      type: record
      fields:
        chromosomes:
          type: File
          inputBinding:
            prefix: --chromosomes
        genomebam:
          type: boolean
          inputBinding:
            prefix: --genomebam
        gtf:
          type: File
          inputBinding:
            prefix: --gtf
  Index:
    type: File
    inputBinding:
      prefix: --index
      position: 1
  InputReads:
    type: File[]
    format: edam:format_1930
    inputBinding:
      position: 200
  PseudoBam:
    type: boolean?
    inputBinding:
      prefix: --pseudobam
  QuantOutfolder:
    type: string
  Seed:
    type: int?
    inputBinding:
      prefix: --seed
  StandardDeviation:
    type: double?
    inputBinding:
      prefix: --sd
  Strand:
    type:
    - 'null'
    - name: forward
      type: record
      fields:
        forward:
          type: boolean
          inputBinding:
            prefix: --fr-stranded
    - name: reverse
      type: record
      fields:
        reverse:
          type: boolean
          inputBinding:
            prefix: --rf-stranded
  isBias:
    type: boolean?
    inputBinding:
      prefix: --bias
  isFusion:
    type: boolean?
    inputBinding:
      prefix: --fusion
  isSingle:
    type: boolean
    inputBinding:
      prefix: --single
      position: 2
  isSingleOverhang:
    type: boolean?
    inputBinding:
      prefix: --single-overhang

outputs:
  kallistoQuantOutDir:
    type: Directory
    outputBinding:
      glob: $(runtime.outdir)/$(inputs.QuantOutfolder)

baseCommand:
- kallisto
- quant
arguments:
- --output-dir
- $(inputs.QuantOutfolder)

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1
  SoftwareRequirement:
    packages:
    - package: kallisto
      specs:
      - https://identifiers.org/biotools/kallisto
      version:
      - 0.51.1

$namespaces:
  edam: https://edamontology.org/
  s: https://schema.org/
$schemas:
- https://edamontology.org/EDAM_1.25.owl
- https://schema.org/version/latest/schemaorg-current-https.rdf
s:citation: https://dx.doi.org/10.1038/nbt.3519
s:codeRepository: https://github.com/pachterlab/kallisto
s:license: https://spdx.org/licenses/BSD-2-Clause
