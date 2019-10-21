Goals: to collect and collaboratively maintain [CWL](https://www.commonwl.org) `CommandLineTool` descriptions of any biology/life-sciences related applications.

[![All Contributors](https://img.shields.io/badge/all_contributors-3-orange.svg?style=flat-square)](#contributors-)
[![Build Status](https://travis-ci.com/common-workflow-library/bio-cwl-tools.svg?branch=release)](https://travis-ci.com/common-workflow-library/bio-cwl-tools)

Non-goals: software packaging or containerization, go to https://biocontainers.pro for that

All CWL tool descriptions are licensed under the Apache 2.0 license.
The underlying tools are under one or more Free and Open Source Software licenses.


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

```yaml
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
```
The first line allows the tool to be run as a single command.
The second specifies the cwl version.
<br/>

**IMPORTANT!** 
The current spec is at v1.1, the IIDSGT group commits tools at **v1.0** due to a limit in our job management system. It is recommended to use the current spec whenever possible.


## Making Files Executable ✴️

Files should be marked as executable before being added 

`
chmod +x tool.cwl
`

## Requirements & Hints Section 🧾

There is a requirements section which handles settings for the runner config. Docker containers should be from biocontainers.pro if possible and placed in the hints section.

```yaml
requirements:
  InlineJavascriptRequirement: {}
```

```yaml
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

If you use schema.org annotations, specify the schema using the RDF version: `$schemas: [ http://schema.org/version/latest/schema.rdf ]`


## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/KerstenBreuer"><img src="https://avatars3.githubusercontent.com/u/28008309?v=4" width="100px;" alt="KerstenBreuer"/><br /><sub><b>KerstenBreuer</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=KerstenBreuer" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/tetron"><img src="https://avatars3.githubusercontent.com/u/1316612?v=4" width="100px;" alt="Peter Amstutz"/><br /><sub><b>Peter Amstutz</b></sub></a><br /><a href="#ideas-tetron" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/michael-kotliar"><img src="https://avatars1.githubusercontent.com/u/19493721?v=4" width="100px;" alt="Michael Kotliar"/><br /><sub><b>Michael Kotliar</b></sub></a><br /><a href="#ideas-michael-kotliar" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=michael-kotliar" title="Code">💻</a></td>
  </tr>
</table>

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
