# gDNAx: Diagnostics for assessing genomic DNA contamination in RNA-seq data

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

install_github("functionalgenomics/gDNAx")
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
If you want to contribute to the development of GenomicScores please open an
[issue](https://github.com/functionalgenomics/gDNAx/issues) to start discussing
your suggestion or, in case of a bugfix or a straightforward feature, directly a
[pull request](https://github.com/functionalgenomics/gDNAx/pulls).
