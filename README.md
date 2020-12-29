<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-12-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END --> 
[![Build Status](https://travis-ci.com/common-workflow-library/bio-cwl-tools.svg?branch=release)](https://travis-ci.com/common-workflow-library/bio-cwl-tools)

Goals: to collect and collaboratively maintain [CWL](https://www.commonwl.org) `CommandLineTool` descriptions of any biology/life-sciences related applications.

Non-goals: software packaging or containerization, go to https://biocontainers.pro for that

All CWL tool descriptions are licensed under the Apache 2.0 license.
The underlying tools are under one or more Free and Open Source Software licenses.

# How to donate your tool descriptions

We welcome pull requests to this repository to add new CWL tool descriptions!

For tools we already have, please compare your definition to the existing
descriptions, maybe you can help improve them.

If you are using a tool (or its subcommand) for a different purpose than the
existing description in this repository, then please submit it as a separate
CWL file next to the existing description(s).

Check out [CONTRIBUTING.md](CONTRIBUTING.md) for a Styleguide and other
contribution tips. Thank you for your contribution!

# How to use these descriptions in your own repository

Here are four different ways that you can use tools from this repository in your own workflows:

1. Make add a [git submodule](https://github.blog/2016-02-01-working-with-submodules/)
   of [this repository](https://github.com/common-workflow-library/bio-cwl-tools) to
   the Git repository of your workflow, typically under the path `tools`.

   [Example of using git submodule with this repo.](https://github.com/arvados/bh20-seq-resource/tree/80cfaba31a99d0c34722312c1b1a69a139477510/workflows)

   Then you can control the exact version of this repository used in your workflow.

2. Copy the entire contents of this repo to your repository. You may have a
   harder time managing updates, but if `git submodule` is uncomfortable
   then this might be a good choice.

   [Example of copying the this repo](https://github.com/mdibl/biocore_analysis/tree/531ae1848cae08c3355175ef3abb048774df866a/biocore_pipelines/cwl_source)

3. Use the in-development [CWL Dependency Manager](https://github.com/common-workflow-language/cwldep).

4. Refer to the tools in this repository by URL, as in this
   [example](https://github.com/pvanheus/lukasa/blob/e0b6c262e1b31521c3a2ea125baf510729cf8950/protein_evidence_mapping.cwl#L29)
   which uses a
   [namespace](https://github.com/pvanheus/lukasa/blob/e0b6c262e1b31521c3a2ea125baf510729cf8950/protein_evidence_mapping.cwl#L89)
   to refer to the  `bio-cwl-tools` repository and then specifies in the individual tool by path within that namespace prefix.
   If this route is used, `cwltool --pack` can created a runnable version of the workflow with all remote references
   resolved.

# Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/KerstenBreuer"><img src="https://avatars3.githubusercontent.com/u/28008309?v=4" width="100px;" alt=""/><br /><sub><b>KerstenBreuer</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=KerstenBreuer" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/tetron"><img src="https://avatars3.githubusercontent.com/u/1316612?v=4" width="100px;" alt=""/><br /><sub><b>Peter Amstutz</b></sub></a><br /><a href="#ideas-tetron" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=tetron" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/michael-kotliar"><img src="https://avatars1.githubusercontent.com/u/19493721?v=4" width="100px;" alt=""/><br /><sub><b>Michael Kotliar</b></sub></a><br /><a href="#ideas-michael-kotliar" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=michael-kotliar" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/pvanheus"><img src="https://avatars0.githubusercontent.com/u/4154788?v=4" width="100px;" alt=""/><br /><sub><b>pvanheus</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=pvanheus" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/medcelerate"><img src="https://avatars3.githubusercontent.com/u/32549017?v=4" width="100px;" alt=""/><br /><sub><b>medcelerate</b></sub></a><br /><a href="#ideas-medcelerate" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/stain"><img src="ihttps://avatars3.githubusercontent.com/u/253413?v=4" width="100px;" alt=""/><br /><sub><b>stain</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=stain" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mb1069"><img src="https://avatars1.githubusercontent.com/u/9156542?v=4" width="100px;" alt=""/><br /><sub><b>Miguel Boland</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=mb1069" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/fjrmoreews"><img src="https://avatars0.githubusercontent.com/u/15047744?v=4" width="100px;" alt=""/><br /><sub><b>fjrmoreews</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=fjrmoreews" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/hellymac"><img src="https://avatars3.githubusercontent.com/u/25847234?v=4" width="100px;" alt=""/><br /><sub><b>cjuigne</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=hellymac" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mberacochea"><img src="https://avatars3.githubusercontent.com/u/1123897?v=4" width="100px;" alt=""/><br /><sub><b>Martín Beracochea</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=mberacochea" title="Code">💻</a></td>
    <td align="center"><a href="https://orcid.org/0000-0002-6130-1021"><img src="https://avatars0.githubusercontent.com/u/1730679?v=4" width="100px;" alt=""/><br /><sub><b>Denis Yuen</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=denis-yuen" title="Documentation">📖</a></td>
    <td align="center"><a href="https://www.linkedin.com/in/sehrish-kanwal-1b80bb42/"><img src="https://avatars3.githubusercontent.com/u/9857259?v=4" width="100px;" alt=""/><br /><sub><b>Sehrish Kanwal</b></sub></a><br /><a href="https://github.com/common-workflow-library/bio-cwl-tools/commits?author=skanwal" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
