#' gDNAx: diagnostics for assessing genomic DNA contamination in RNA-seq data
#' 
#' The gDNAx package provides diagnostics for assessing genomic DNA
#' contamination in RNA-seq data, as well as plots representing these
#' diagnostics. Moreover, the package can be used to get an insight into
#' the strand library protocol used and, in case of strand-specific libraries,
#' the strandedness of the data. Furthermore, it provides functionality to
#' filter out reads of potential gDNA origin.
#' 
#' The main functions are:
#'
#' \itemize{
#'     \item \code{\link{gDNAdx}()} - calculate diagnostics for assessing the presence of genomic DNA in RNA-seq data over a subset of the alignments in the input BAM files.
#'     \item \code{\link{getDx}()} and \code{\link{plot}()} - get and plot statistics on genomic DNA contamination levels, respectively.
#'     \item \code{\link{strandedness}()} - obtain estimates of strandedness in RNA-seq data samples based on the proportion of reads aligning to the same or opposite strand as transcripts in the annotations.
#'     \item \code{\link{classifyStrandMode}()} - classify the output of \code{\link{strandedness}()} into strand modes for each BAM file.
#'     \item \code{\link{filterBAMtxFlag}} and \code{\link{filterBAMtx}} - filter alignments in a BAM file using criteria based on a transcriptome annotation.
#' }
#' 
#' For detailed information on usage, see the package vignette, by typing
#' \code{vignette("gDNAx")}.
#' 
#' All questions and bug reports should be posted to the Bioconductor Support 
#' Site:
#' 
#' \url{https://support.bioconductor.org}
#'
#' The code of the development version of the package is available at the
#' GitHub repository:
#'
#' \url{https://github.com/functionalgenomics/gDNAx}
#' 
#' @name gDNAx-package
#' @aliases gDNAx-package
#' @aliases gDNAx
#' @keywords package
"_PACKAGE"

NULL
