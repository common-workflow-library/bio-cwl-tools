#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

doc: |
  Performs strand cross-correlation analysis of aligned reads and estimates the fragment size.

requirements:
  - class: InlineJavascriptRequirement
hints:
  DockerRequirement:
    dockerPull: kerstenbreuer/phantompeakqualtools:1.2
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000

# Please note: please adjust the path to run_spp.R if not using containers:
baseCommand: ["Rscript", "--verbose", "--max-ppsize=500000", "/usr/bin/phantompeakqualtools-1.2/run_spp.R"]
arguments:
  - valueFrom: $(runtime.tmpdir)
    prefix: -tmpdir=
    separate: false
    position: 10
  - valueFrom: $(runtime.outdir)
    prefix: -odir=
    separate: false
    position: 10
  - valueFrom: $(inputs.bam.nameroot).crosscor.pdf
    prefix: -savp=
    separate: false
    position: 100
  - valueFrom: $(inputs.bam.nameroot).spp.out
    prefix: -out=
    separate: false
    position: 100
stderr: $(inputs.bam.nameroot).phantompeakqualtools_stderr
stdout: $(inputs.bam.nameroot).phantompeakqualtools_stdout # Naming important for multiqcs
 
inputs:
  bam:
    type: File
    inputBinding:
      prefix: -c=
      separate: false
      position: 10
 
outputs:
  qc_crosscorr_summary:
    type: File?
    outputBinding:
      glob:  "*.spp.out"
  qc_crosscorr_plot:
    type: File?
    outputBinding:
      glob:  "*.pdf"
  qc_crosscorr_fragment_size:
    type: int?
    outputBinding:
      glob:  "*.spp.out"
      loadContents: true
      outputEval: $(parseInt(self[0].contents.split('\t')[2]))
  qc_phantompeakqualtools_stderr:
    type: stderr
  qc_phantompeakqualtools_stdout:
    type: stdout
