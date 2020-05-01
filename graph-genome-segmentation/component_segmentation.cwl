#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  ResourceRequirement:
    coresMin: 8
  DockerRequirement:
    dockerImageId: graph-genome-component_segmentation
    dockerFile: |
        FROM debian:stable-slim
        RUN apt-get update && apt-get install -y git python3-pip && rm -rf /var/lib/apt/lists/*
        WORKDIR /src
        RUN git clone --depth 1 https://github.com/graph-genome/component_segmentation && rm -Rf component_segmentation/data
        WORKDIR /src/component_segmentation
        RUN pip3 install --no-cache-dir -r requirements.txt
        RUN chmod a+x segmentation.py
        ENV PATH="/src/component_segmentation:${PATH}"


inputs:
  bins:
    type: File
    inputBinding:
      prefix: --json-file=
      separate: false

  cells_per_file:
    type: int?
    label: Number of cells per file
    doc: |
      #number bins per file = #cells / #paths
      Tip: Adjust this number to get chunk files output close to 2MB.
    inputBinding:
      prefix: --cells-per-file=
      separate: false

arguments:
  - --out-folder=$(runtime.outdir)
  - --parallel-cores=$(runtime.cores)

baseCommand: segmentation.py

outputs:
  colinear_components:
    type: File[]
    outputBinding:
      glob: "*/*.schematic.json"
