# Styleguide

## Naming Tools 📛

Tools should follow the convention of being prefixed by the parent tool name and Camelcase like so i.e.
`
BWA-Mem.cwl`
or `
BWA-Index.cwl
`

## Tool Feature Requirements 🆕

The first 3 lines of tool wrappers should be as follows. Our CI/CD system checks for these so make sure to include them so they can be merged into the repo.

``` cwl
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
```
The first line allows the tool to be run as a single command.
The second specifies the cwl version.
<br/>

**IMPORTANT!** 
We are using `cwlVersion` `v1.0` unless a `v1.1` feature is needed.

## Making Files Executable ✴️

Files should be marked as executable before being added 

`
chmod +x tool.cwl
`

## Requirements & Hints Section 🧾

There is a requirements section which handles settings for the runner config. Docker containers should be from biocontainers.pro if possible and placed in the hints section.

``` cwl
requirements:
  InlineJavascriptRequirement: {}
```

``` cwl
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
```

## Validation ✅

Tools need to be free of warning when running with

`
cwltool --validate
`

## Adding Tools To The Repository ➕

Please add tools via pull requests to the repository. Our CI/CD runs validation against the tools and will soon support doing unit tests on the individual tools.

## Descriptions 📃

Tool descriptions should be motivated by a real world use of this tool in a workflow.
The description should focus on a single way of using the tool.
Signs that a tool description is including too much: lots of javascript; complicated data structures; every single flag is listed.

## Schema Description

If you use schema.org annotations, specify the schema using the RDF version:
`$schemas: [ http://schema.org/version/latest/schema.rdf ]` unless items from
outside the core schema.org vocabulary are needed. In that case use
`$schemas: [ https://schema.org/version/latest/all-layers.rdf ]`.

However, don't use `s:mainEntity`, put that information under `hints` as a `SoftwareRequirement`.

## File Formats

If your tool has well defined input or output files, we recommend the addition of file formats using ontologies such as EDAM. CWL executors like cwltool can do some basic reasoning using this information and can warn about obvious mismatches.  

``` cwl
input_sequences:
    type: 
      - File
      - File[]
    label: "Input sequence files"
    format:
      - edam:format_1929  # FASTA
      - edam:format_1930  # FASTQ
```
