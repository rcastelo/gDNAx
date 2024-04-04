#' Remove gDNA contamination from RNA-seq data
#'
#' Remove gDNA contamination from RNA-seq data by filtering read alignments in
#' BAM files that putatively have a gDNA origin. This is currently a wrapper
#' with convenient default values for the function \code{\link{filterBAMtx}()},
#' please use that function if you need greater control on how to filter RNA-seq
#' alignments.
#'
#' @param x \linkS4class{gDNAx} object obtained with the function
#' \code{\link{gDNAdx}()}.
#'
#' @param path Directory where to write the filtered BAM files.
#'
#' @param sbparam Either \code{NULL} (default) or a \linkS4class{ScanBamParam}
#' object. The \code{NULL} value implies that internally a
#' \linkS4class{ScanBamParam} object is built with the following flags:
#' \code{isUnmappedQuery=FALSE}, \code{isProperPair=!singleEnd(x)},
#' \code{isSecondaryAlignment=FALSE}, \code{isNotPassingQualityControls=FALSE},
#' \code{isDuplicate=FALSE}.
#'
#' @param yieldSize (Default 1e6) Number of records in the input BAM file to
#' yield each time the file is read. The lower the value, the smaller memory
#' consumption, but in the case of large BAM files, values below 1e6 records
#' may decrease the overall performance.
#'
#' @param verbose (Default TRUE) Logical value indicating if progress should be
#' reported through the execution of the code.
#'
#' @param BPPARAM An object of a \linkS4class{BiocParallelParam} subclass
#' to configure the parallel execution of the code. By default, a
#' \linkS4class{SerialParam} object is used, which does not use any
#' parallelization, with the flag \code{progress=TRUE} to show progress
#' through the calculations.
#'
#' @return A \code{data.frame} object with the number of filtered read alignments
#' tallied by their origin.
#'
#' @examples
#'
#' library(gDNAinRNAseqData)
#'  
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' ## fetch sample BAM files
#' bamfiles <- LiYu22subsetBAMfiles()
#' bamfiles <- bamfiles[c(1,7)] # using a subset of samples
#'
#' ## diagnose gDNA contamination
#' gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#'
#' ## remove gDNA contamination
#' dir <- tempdir()
#' fstats <- gDNAtx(gdnax, path=dir)
#' fstats
#' list.files(dir, pattern="*.bam$")
#'
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom BiocParallel SerialParam
#' @export
#' @rdname gDNAtx
gDNAtx <- function(x, path=".", sbparam=NULL, yieldSize=1000000L, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
    if (!is(x, "gDNAx"))
        stop("'x' should be a 'gDNAx' object.")
    if (!file.exists(path))
        stop(sprintf("path %s does not exist.", path))

    yieldSize <- .checkYieldSize(yieldSize)

    sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                           isProperPair=!singleEnd(x),
                           isSecondaryAlignment=FALSE,
                           isNotPassingQualityControls=FALSE,
                           isDuplicate=FALSE)
    if (is.null(sbparam))
        sbparam <- ScanBamParam(flag=sbflags)
    else {
        if (!is(sbparam, "ScanBamParam"))
            stop("'sbparam' should be either 'NULL' or a 'ScanBamParam' object.")
    }

    fbf <- filterBAMtxFlag(isSpliceCompatibleJunction=TRUE,
                           isSpliceCompatibleExonic=TRUE,
                           isInStrandedWindow=TRUE)

    out <- filterBAMtx(x, path=path, txflag=fbf, param=sbparam,
                       yieldSize=yieldSize, verbose=verbose,
                       BPPARAM=BPPARAM)
    out
}
