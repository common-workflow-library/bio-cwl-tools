Goals: to collect and collaboratively maintain [CWL](https://www.commonwl.org) `CommandLineTool` descriptions of any biology/life-sciences related applications.
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-9-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END --> 
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
We are using `cwlVersion` `v1.0` unless a `v1.1` feature is needed.

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

However, don't use `s:mainEntity`, put that information under `hints` as a `SoftwareRequirement`

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/KerstenBreuer"><img src="https://avatars3.githubusercontent.com/u/28008309?v=4" width="100px;" alt=""/><br /><sub><b>KerstenBreuer</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=KerstenBreuer" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/tetron"><img src="https://avatars3.githubusercontent.com/u/1316612?v=4" width="100px;" alt=""/><br /><sub><b>Peter Amstutz</b></sub></a><br /><a href="#ideas-tetron" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/michael-kotliar"><img src="https://avatars1.githubusercontent.com/u/19493721?v=4" width="100px;" alt=""/><br /><sub><b>Michael Kotliar</b></sub></a><br /><a href="#ideas-michael-kotliar" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=michael-kotliar" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/pvanheus"><img src="https://avatars0.githubusercontent.com/u/4154788?v=4" width="100px;" alt=""/><br /><sub><b>pvanheus</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=pvanheus" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/medcelerate"><img src="https://avatars3.githubusercontent.com/u/32549017?v=4" width="100px;" alt=""/><br /><sub><b>medcelerate</b></sub></a><br /><a href="#ideas-medcelerate" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/stain"><img src="ihttps://avatars3.githubusercontent.com/u/253413?v=4" width="100px;" alt=""/><br /><sub><b>stain</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=stain" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mb1069"><img src="https://avatars1.githubusercontent.com/u/9156542?v=4" width="100px;" alt=""/><br /><sub><b>Miguel Boland</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=mb1069" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/fjrmoreews"><img src="https://avatars0.githubusercontent.com/u/15047744?v=4" width="100px;" alt=""/><br /><sub><b>fjrmoreews</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=fjrmoreews" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/hellymac"><img src="https://avatars3.githubusercontent.com/u/25847234?v=4" width="100px;" alt=""/><br /><sub><b>cjuigne</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=hellymac" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
