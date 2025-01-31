# gDNAx: Diagnostics for assessing genomic DNA contamination in RNA-seq data

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/gDNAx.svg)](https://bioconductor.org/packages/release/bioc/html/gDNAx.html "How long has been gDNAx in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/release/gDNAx.svg)](https://bioconductor.org/packages/stats/bioc/gDNAx/ "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months.")
[![Support posts](https://bioconductor.org/shields/posts/gDNAx.svg)](https://support.bioconductor.org/t/gDNAx/ "Support site activity on gDNAx, last 6 months: answered posts/total posts.")
[![R-CMD-check-bioc](https://github.com/rcastelo/gDNAx/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/rcastelo/gDNAx/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/rcastelo/gDNAx/coverage.svg?branch=devel)](https://app.codecov.io/github/rcastelo/gDNAx?branch=devel)
<img align="right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/gDNAx/gDNAx.png" height="200"/>

**Current Bioconductor build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/gDNAx.svg)](https://bioconductor.org/packages/release/bioc/html/gDNAx.html#archives "Whether gDNAx release is available on all platforms")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/gDNAx.svg)](https://bioconductor.org/packages/release/bioc/html/gDNAx.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/gDNAx.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/gDNAx "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/gDNAx.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/gDNAx/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/gDNAx.svg)](https://bioconductor.org/packages/devel/bioc/html/gDNAx.html#archives "Whether gDNAx devel is available on all platforms")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/gDNAx.svg)](https://bioconductor.org/packages/devel/bioc/html/gDNAx.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/gDNAx.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/gDNAx "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/gDNAx.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/gDNAx/ "Bioconductor devel build")

The `gDNAx` package calculates and plots different diagnostics related to gDNA
contamination levels, such as strandedness and proportion of intergenic reads,
intronic reads, splice compatible exonic reads and splice compatible junction
reads. Moreover, the package is able to identify the strand library type in
stranded RNA-seq data.

## Installation

You can install the `gDNAx` package from this GitHub repo using
the [remotes](https://cran.r-project.org/package=remotes) package, as follows:

```
library(remotes)

install_github("rcastelo/gDNAx")
```

Provided that you have installed first all its Bioconductor dependencies;
see the `DESCRIPTION` file. The vignette contains an example on how to use it.

## Questions, bug reports and issues

For bug reports and issues regarding this __development__ version of **gDNAx**
please use the GitHub issues [tab](https://github.com/functionalgenomics/gDNAx/issues)
at the top-left of this page.

## Contributing

Contributions to the software codebase of gDNAx are welcome as long as contributors
abide to the terms of the 
[Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct).
If you want to contribute to the development of gDNAx please open an
[issue](https://github.com/functionalgenomics/gDNAx/issues) to start discussing
your suggestion or, in case of a bugfix or a straightforward feature, directly a
[pull request](https://github.com/functionalgenomics/gDNAx/pulls).

## Funding

This software project is supported in part by the Spanish
[Ministry of Science, Innovation and Universities](https://www.ciencia.gob.es).
