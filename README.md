Goals: to collect and collaboratively maintain [CWL](https://www.commonwl.org) `CommandLineTool` descriptions of any biology/life-sciences related applications.

[![All Contributors](https://img.shields.io/badge/all_contributors-3-orange.svg?style=flat-square)](#contributors)

Non-goals: software packaging or containerization, go to https://biocontainers.pro for that

All CWL tool descriptions are licensed under the Apache 2.0 license.
The underlying tools are under one or more Free and Open Source Software licenses.


# Styleguide

1. First line must be `#!/usr/bin/env cwl-wrapper`
1. Second line is `cwlVersion: v1.0`
1. Third line is `class: CommandLineTool`
1. `*.cwl` files must be marked executable (`chmod a+x *.cwl`)
1. All tool descriptions must have a software container. Use a container from biocontainers.pro if available
1. If you use schema.org annotations, specify the schema using the RDF version: `$schemas: [ http://schema.org/version/latest/schema.rdf ]`

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
