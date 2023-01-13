#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement:
    expressionLib:
    - var get_label = function(i) {
          var rootname = inputs.molecule_info_h5[i].basename.split('.').slice(0,-1).join('.');
          rootname = (rootname=="")?inputs.molecule_info_h5[i].basename:rootname;
          return inputs.gem_well_labels?inputs.gem_well_labels[i].replace(/,/g, "_"):rootname;
      };
  InitialWorkDirRequirement:
    listing: |
      ${
        var entry = "library_id,molecule_h5\n"
        for (var i=0; i < inputs.molecule_info_h5.length; i++){
          entry += get_label(i) + "," + inputs.molecule_info_h5[i].path + "\n"
        }
        return [{
          "entry": entry,
          "entryname": "metadata.csv"
        }];
      }
hints:
  DockerRequirement:
    dockerPull: cumulusprod/cellranger:4.0.0

inputs:  
  molecule_info_h5:
    type: File[]
    format: edam:format_3590  # HDF5
    doc: |
      Array of molecule-level information files in HDF5 format.
      Outputs from "cellranger count" command
  
  gem_well_labels:
    type:
    - "null"
    - string[]
    doc: |
      Array of GEM well identifiers to be used for labeling purposes only.
      If not provided use rootnames of files from the molecule_info_h5 input

  normalization_mode:
    type:
    - "null"
    - type: enum
      name: "normalization"
      symbols: ["none", "mapped"]
    inputBinding:
      position: 5
      prefix: "--normalize"
    doc: |
      Library depth normalization mode: mapped, none.
      Default: mapped

  threads:
    type: int?
    inputBinding:
      position: 6
      prefix: "--localcores"
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available

  memory_limit:
    type: int?
    inputBinding:
      position: 7
      prefix: "--localmem"
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available

  virt_memory_limit:
    type: int?
    inputBinding:
      position: 8
      prefix: "--localvmem"
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available


outputs:

  web_summary_report:
    type: File
    format: iana:text/html
    outputBinding:
      glob: "aggregated/outs/web_summary.html"
    doc: |
      Aggregated run summary metrics and charts in HTML format

  metrics_summary_report_json:
    type: File
    format: iana:application/json
    outputBinding:
      glob: "aggregated/outs/summary.json"
    doc: |
      Aggregated run summary metrics in JSON format
  
  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/analysis"
    doc: |
      Folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression for aggregated results

  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/filtered_feature_bc_matrix"
    doc: |
      Folder with aggregated filtered feature-barcode matrices containing only cellular barcodes in MEX format

  filtered_feature_bc_matrix_h5:
    type: File
    format: edam:format_3590  # HDF5
    outputBinding:
      glob: "aggregated/outs/filtered_feature_bc_matrix.h5"
    doc: |
      Aggregated filtered feature-barcode matrices containing only cellular barcodes in HDF5 format

  raw_feature_bc_matrices_folder:
    type: Directory
    outputBinding:
      glob: "aggregated/outs/raw_feature_bc_matrix"
    doc: |
      Folder with aggregated unfiltered feature-barcode matrices containing all barcodes in MEX format

  raw_feature_bc_matrices_h5:
    type: File
    format: edam:format_3590  # HDF5
    outputBinding:
      glob: "aggregated/outs/raw_feature_bc_matrix.h5"
    doc: |
      Aggregated unfiltered feature-barcode matrices containing all barcodes in HDF5 format

  aggregation_metadata:
    type: File
    format: iana:text/csv
    outputBinding:
      glob: "aggregated/outs/aggregation.csv"
    doc: |
      Copy of the input aggregation CSV file

  loupe_browser_track:
    type: File
    outputBinding:
      glob: "aggregated/outs/cloupe.cloupe"
    doc: |
      Loupe Browser visualization and analysis file for aggregated results

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger", "aggr", "--disable-ui", "--id", "aggregated", "--csv", "metadata.csv"]


stdout: cellranger_aggr_stdout.log
stderr: cellranger_aggr_stderr.log


$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/
  iana: https://www.iana.org/assignments/media-types/

$schemas:
- https://schema.org/version/latest/schemaorg-current-https.rdf

label: "Cellranger aggr - aggregates data from multiple Cellranger runs"
s:alternateName: "Cellranger aggr takes a list of cellranger count output files and produces a single feature-barcode matrix containing all the data"

s:license: http://www.apache.org/licenses/LICENSE-2.0

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Tool calls "cellranger aggr" command to combine output files from "cellranger count"
  (the molecule_info.h5 file from each run) into a single feature-barcode matrix containing
  all the data. When combining multiple GEM wells, the barcode sequences for each channel
  are distinguished by a GEM well suffix appended to the barcode sequence. Each GEM well is
  a physically distinct set of GEM partitions, but draws barcode sequences randomly from the
  pool of valid barcodes, known as the barcode whitelist. To keep the barcodes unique when
  aggregating multiple libraries, we append a small integer identifying the GEM well to the
  barcode nucleotide sequence, and use that nucleotide sequence plus ID as the unique identifier
  in the feature-barcode matrix. For example, AGACCATTGAGACTTA-1 and AGACCATTGAGACTTA-2 are
  distinct cell barcodes from different GEM wells, despite having the same barcode nucleotide
  sequence. This number, which tells us which GEM well this barcode sequence came from, is
  called the GEM well suffix. The numbering of the GEM wells will reflect the order that the
  GEM wells were provided in the "molecule_info_h5" and "gem_well_labels" inputs.

  When combining data from multiple GEM wells, the "cellranger aggr" pipeline automatically
  equalizes the average read depth per cell between groups before merging. This approach avoids
  artifacts that may be introduced due to differences in sequencing depth. It is possible to turn
  off normalization or change the way normalization is done through the "normalization_mode"
  input. The "none" value may be appropriate if you want to maximize sensitivity and plan to deal
  with depth normalization in a downstream step.

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id - hardcoded to `aggregated` as we want to return the content of the
         outputs folder as separate outputs

  Skipped parameters:
  --nosecondary
  --dry
  --noexit
  --nopreflight
  --description
  --jobmode
  --mempercore
  --maxjobs
  --jobinterval
  --overrides
  --uiport

  Not supported features:
  - Batch correction caused by different versions of the Single Cell Gene Expression chemistry is
    not supported as the generated metadata file doesn't include "batch" field.
