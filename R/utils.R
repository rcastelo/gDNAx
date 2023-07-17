## private function .checkBamFileListArgs()
## adapted from GenomicAlignments/R/summarizeOverlaps-methods.R

#' @importFrom methods is
#' @importFrom Rsamtools BamFileList asMates asMates<-
.checkBamFileListArgs <- function(bfl, singleEnd, fragments, yieldSize) {
    if (missing(bfl) || length(bfl) == 0 ||
        !class(bfl) %in% c("character", "BamFileList"))
        stop("argument 'bfl' should be either a string character vector",
                "of BAM file names or a 'BamFileList' object")
    
    if (is.character(bfl)) {
        mask <- vapply(bfl, FUN = file.exists, FUN.VALUE = logical(1))
        if (any(!mask)) {
            whmiss <- paste(paste("  ", bfl[!mask]), collapse="\n")
            stop(sprintf("The following input BAM files cannot be found:\n%s",
                        whmiss))
        }
        if (any(duplicated(bfl))) {
            whdupl <- paste(paste("  ", bfl[duplicated(bfl)]), collapse="\n")
            stop(sprintf("The following input BAM files are duplicated:\n%s",
                        whdupl))
        }
    }
    
    if (!is(bfl, "BamFileList"))
        bfl <- BamFileList(bfl, asMates=!singleEnd, yieldSize=yieldSize)
    
    if (singleEnd) {
        if (all(isTRUE(asMates(bfl))))
            stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
        if (fragments)
            stop("when 'fragments=TRUE', 'singleEnd' should be FALSE")
    } else
        asMates(bfl) <- TRUE
    
    bfl
}

## copied from GenomicAlignments:::.normarg_strandMode_replace_value()
#' @importFrom S4Vectors isSingleNumber
.checkStrandMode <- function(value) {
    if (!is.na(value)) {
        if (!isSingleNumber(value))
            stop("invalid strand mode (must be 0, 1, or 2)")
        if (!is.integer(value))
            value <- as.integer(value)
        if (!(value %in% 0:2))
            stop("invalid strand mode (must be 0, 1, or 2)")
    } else {
        value <- as.integer(value)
    }

    value
}

## private function .checkYieldSize()
.checkYieldSize <- function(value) {
    if (is.na(value) || !is.numeric(value) || value < 10000) {
        errmsg <- paste("the parameter 'yieldSize' should be a number larger",
                        "or equal than 10000")
        stop(errmsg)
    }
    value <- as.integer(value)
    value
}

## private function .checkPairedEnd()
## checks whether BAM files are all single end
## or paired end, and whether that matches the
## input argument 'singleEnd'

#' @importFrom Rsamtools yieldSize yieldSize<- testPairedEndBam
#' @importFrom S4Vectors unname
.allPairedEnd <- function(bfl, singleEnd) {
    testpe <- function(bf) {
        yieldSize <- yieldSize(bf)
        yieldSize(bf) <- 10000 ## use a small yieldSize for speed
        on.exit(yieldSize(bf) <- yieldSize)
        on.exit(close(bf), add=TRUE)
        if (!isOpen(bf))
            bf <- open(bf)
        testPairedEndBam(bf)
    }
    peflag <- unname(vapply(bfl, testpe, FUN.VALUE = logical(1L)))
    if (singleEnd && !all(!peflag))
        stop("Some BAM files are paired end, but 'singleEnd=TRUE'.")

    if (!singleEnd && !all(peflag))
        stop("Some BAM files are single end, but 'singleEnd=FALSE'.")

    peflag[1]
}
## private function .minFrgLen()
## returns the read length for single-end and
## twice the read length for paired-end

#' @importFrom Rsamtools yieldSize yieldSize<- scanBam
#' @importFrom S4Vectors unname
.minFrgLen <- function(bfl, singleEnd) {
    stopifnot(length(bfl) > 0) ## QC
    rlen <- rep(0L, length(bfl))
    peflag <- .allPairedEnd(bfl, singleEnd)
    qw <- function(bf) {
        yieldSize <- yieldSize(bf)
        yieldSize(bf) <- 10000 ## use a small yieldSize for speed
        on.exit(yieldSize(bf) <- yieldSize)
        on.exit(close(bf), add=TRUE)
        if (!isOpen(bf))
            bf <- open(bf)
        w <- scanBam(bf, param=ScanBamParam(what="qwidth"))[[1]]$qwidth
        tab <- table(w)
        w <- as.integer(names(tab[which.max(tab)]))
        w
    }
    allw <- unname(vapply(bfl, qw, FUN.VALUE = integer(1L)))
    if (length(unique(allw)) > 1)
        warning(sprintf("BAM files have different read lengths, using the minimum (%d).", min(allw)))

    if (peflag)
        allw <- allw * 2

    min(allw)
}

## private function .ppprintnames() for pretty-printing
## a vector of character strings
.pprintnames <- function(x) {
    y <- x
    if (length(x) > 2)
        y <- c(y[1], "...", y[length(y)])
    y <- paste(y, collapse=", ")
    y
}

## private function .loadAnnotationPackageObject()
#' @importFrom utils installed.packages
.loadAnnotationPackageObject <- function(pkgName, argName, pkgType,
                                         verbose=TRUE) {

    callobj <- match.call()
    annotObj <- NULL

    if (is.character(pkgName)) {
        if (!pkgName %in% installed.packages(noCache=TRUE)[, "Package"])
            stop(sprintf("please install the Bioconductor package %s.", 
                         pkgName))
        if (!.isPkgLoaded(pkgName)) {
            if (verbose)
                message("Loading ", pkgType, " annotation package ", pkgName)
            # loaded <- suppressPackageStartupMessages(require(pkgName,
            #                                             character.only=TRUE))
            loaded <- suppressPackageStartupMessages(requireNamespace(pkgName,
                                                        character.only=TRUE))
            if (!loaded)
                stop(sprintf("package %s could not be loaded.", pkgName))
        }
        tryCatch({
            # annotObj <- get(pkgName)
            attachNamespace(pkgName)
            annotObj <- get(pkgName)
        }, error=function(err) {
            fmtstr <- paste("The annotation package %s should automatically load",
                            "an %s object with the same name as the package.")
            errmsg <- sprintf(fmtstr, pkgName, pkgType)
            stop(errmsg)
        })
    } else if (class(pkgName) != pkgType) {
        fmtstr <- paste("argument '%s' should either contain the name of an",
                        "'%s' annotation package or be an '%s' annotation",
                        "object itself.")
        errmsg <- sprintf(fmtstr, argName, pkgType, pkgType)
        stop(errmsg)
    } else
        annotObj <- pkgName

    if (!is(annotObj, pkgType))
        stop(sprintf("The object loaded with name %s is not an '%s' object.",
                    ifelse(is.character(pkgName), pkgName,
                            gettext(callobj)[2])), pkgType)

    annotObj
}

## private function .isPkgLoaded()
.isPkgLoaded <- function(name) {
    (paste("package:", name, sep="") %in% search())
}

## private function .matchSeqinfo()

#' @importFrom GenomeInfoDb seqlengths keepSeqlevels seqlevelsStyle
#' @importFrom GenomeInfoDb seqinfo seqinfo<- seqlevelsStyle<- seqlevels
.matchSeqinfo <- function(gal, tx, verbose=TRUE) {
    stopifnot("GAlignments" %in% class(gal) ||
                "GAlignmentPairs" %in% class(gal) ||
                "GAlignmentsList" %in% class(gal) ||
                "TxDb" %in% class(tx)) ## QC

    seqlevelsStyle(gal) <- seqlevelsStyle(tx)[1]
    slengal <- seqlengths(gal)
    slentx <- seqlengths(tx)
    commonchr <- intersect(names(slengal), names(slentx))
    slengal <- slengal[commonchr]
    slentx <- slentx[commonchr]
    if (any(slengal != slentx)) {
        if (sum(slengal != slentx) == 1 && verbose) {
            difflonechr <- paste("Chromosome %s has different lengths",
                                "between the input BAM and the TxDb",
                                "annotation package. This chromosome will",
                                "be discarded from further analysis",
                                sep=" ")
            whdifflonechr <- paste(commonchr[which(slengal != slentx)],
                                    collapse=", ")
            message(sprintf(difflonechr, whdifflonechr))
        } else if (verbose) {
            difflmultichr <- paste("Chromosomes %s have different lengths",
                                    "between the input BAM and the TxDb",
                                    "annotation package. These chromosomes",
                                    "will be discarded from further analysis",
                                    sep=" ")
            whdifflmultichr <- paste(commonchr[which(slengal != slentx)],
                                    collapse=", ")
            message(sprintf(difflmultichr, whdifflmultichr))
        }
        if (sum(slengal == slentx) == 0)
            stop("None of the chromosomes in the input BAM file has the ",
                    "same length as the chromosomes in the input TxDb ",
                    "annotation package.")
        gal <- keepSeqlevels(gal, commonchr[slengal == slentx],
                            pruning.mode="coarse")
        commonchr <- commonchr[slengal == slentx]
    }

    ## set the seqinfo information to the one of the TxDb annotations
    mt <- match(commonchr, seqlevels(gal))
    seqinfo(gal, new2old=mt, pruning.mode="coarse") <- seqinfo(tx)[commonchr]

    gal
}

## private function .getReadFunction()
## borrowed from GenomicAlignments/R/summarizeOverlaps-methods.R
#' @importFrom GenomicAlignments readGAlignmentsList readGAlignments 
#' @importFrom GenomicAlignments readGAlignmentPairs
.getReadFunction <- function(singleEnd, fragments) {
    if (singleEnd) {
        FUN <- readGAlignments
    } else {
        if (fragments)
            FUN <- readGAlignmentsList
        else
            FUN <- readGAlignmentPairs
    }
    
    FUN
}

## private function .appendHits()
## appends the second Hits object to the end of the first one
## assuming they have identical right nodes
#' @importFrom S4Vectors nLnode nRnode isSorted from to Hits
.appendHits <- function(hits1, hits2) {
    stopifnot(nRnode(hits1) == nRnode(hits2))
    stopifnot(isSorted(from(hits1)) == isSorted(from(hits2)))
    hits <- c(Hits(from=from(hits1), to=to(hits1),
                    nLnode=nLnode(hits1)+nLnode(hits2),
                    nRnode=nRnode(hits1), sort.by.query=isSorted(from(hits1))),
              Hits(from=from(hits2)+nLnode(hits1), to=to(hits2),
                    nLnode=nLnode(hits1)+nLnode(hits2),
                    nRnode=nRnode(hits2), sort.by.query=isSorted(from(hits2))))
    hits
}
