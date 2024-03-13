## private function .checkBamFileListArgs()

#' @importFrom methods is
#' @importFrom Rsamtools BamFileList
#' @importFrom BiocGenerics path
.checkBamFiles <- function(bfl) {
    if (missing(bfl) || length(bfl) == 0 ||
        !class(bfl) %in% c("character", "BamFileList"))
        stop("argument 'bfl' should be either a string character vector",
             "of BAM file names or a 'BamFileList' object")
    
    fnames <- bfl
    if (is(bfl, "BamFileList"))
      fnames <- path(bfl)

    mask <- vapply(fnames, FUN=file.exists, FUN.VALUE=logical(1))
    if (any(!mask)) {
        whmiss <- paste(paste("  ", bfl[!mask]), collapse="\n")
        stop(sprintf("The following input BAM files cannot be found:\n%s",
                     whmiss))
    }
    if (any(duplicated(fnames))) {
        whdupl <- paste(paste("  ", bfl[duplicated(bfl)]), collapse="\n")
        stop(sprintf("The following input BAM files are duplicated:\n%s",
                     whdupl))
    }
    fsize <- vapply(fnames, FUN=file.size, FUN.VALUE=numeric(1))
    mask <- fsize > 0
    if (any(!mask)) {
        whzero <- paste(paste("  ", bfl[!mask]), collapse="\n")
        stop(sprintf("The following input BAM files have 0 bytes:\n%s",
                     whzero))
    }

    bfl
}
    
## private function .checkBamFileListArgs() it assumes that
## .checkBamFiles() and .checkYieldSize() have been called before

#' @importFrom Rsamtools BamFileList asMates asMates<-
.checkBamFileListArgs <- function(bfl, singleEnd, yieldSize) {
    if (is.character(bfl))
        bfl <- BamFileList(bfl, asMates=!singleEnd, yieldSize=yieldSize)
    
    if (singleEnd) {
        if (all(isTRUE(asMates(bfl))))
            stop("cannot specify both 'singleEnd=TRUE' and 'asMates=TRUE'")
    } else
        asMates(bfl) <- TRUE
    
    bfl
}

## private function .checkStrandMode()

#' @importFrom S4Vectors isSingleNumber
.checkStrandMode <- function(strandMode) {
    if (!is.na(strandMode)) {
        if (!isSingleNumber(strandMode))
            stop("invalid strand mode (must be 0, 1, or 2)")
        if (!is.integer(strandMode))
            strandMode <- as.integer(strandMode)
        if (!(strandMode %in% 0:2))
            stop("invalid strand mode (must be 0, 1, or 2)")
    } else {
        strandMode <- as.integer(strandMode)
    }

    strandMode
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
## checks whether BAM files have a single- or paired-end layout,
## whether they all have the same layout, and if 'singleEnd' is
## not missing, it also checks whether its value matches the
## layout found in the BAM files. it assumes that .checkBamFiles()
## has been called before.

#' @importFrom Rsamtools yieldSize yieldSize<- testPairedEndBam
#' @importFrom BiocParallel bpnworkers
#' @importFrom S4Vectors unname
.checkPairedEnd <- function(bfl, singleEnd, BPPARAM=SerialParam()) {
    if (is.character(bfl))
      bfl <- BamFileList(bfl)
    testpe <- function(bf) {
        yieldSize <- yieldSize(bf)
        yieldSize(bf) <- 10000 ## use a small yieldSize for speed
        on.exit(yieldSize(bf) <- yieldSize)
        on.exit(close(bf), add=TRUE)
        if (!isOpen(bf))
            bf <- open(bf)
        suppressMessages(testPairedEndBam(bf))
    }
    peflag <- FALSE
    if (length(bfl) > 1 && bpnworkers(BPPARAM) > 1) {
      bpparsilent <- BPPARAM
      bpprogressbar(bpparsilent) <- FALSE
      peflag <- unlist(bplapply(bfl, testpe, BPPARAM=bpparsilent),
                       use.names=FALSE)
    } else
      peflag <- unname(vapply(bfl, testpe, FUN.VALUE = logical(1L)))

    if (!all(peflag[1] == peflag))
      stop("Some BAM files are single-end and some are paired-end.")
    if (!missing(singleEnd)) {
      if (singleEnd && any(peflag))
          stop("Some BAM files are paired-end, but 'singleEnd=TRUE'.")

      if (!singleEnd && any(!peflag))
          stop("Some BAM files are single-end, but 'singleEnd=FALSE'.")
    }

    peflag[1]
}

## private function readLengths()
## figure out read length from each BAM file

#' @importFrom Rsamtools yieldSize yieldSize<- scanBamFlag scanBam
#' @importFrom S4Vectors unname
#' @importFrom cli cli_alert_info
.readLengths <- function(bfl, singleEnd, verbose, BPPARAM=SerialParam()) {
    stopifnot(length(bfl) > 0) ## QC
    qw <- function(bf, singleEnd) { ## figure out most freq. query width in first 10K aln.
        yieldSize <- yieldSize(bf)
        yieldSize(bf) <- 10000 ## use a small yieldSize for speed
        on.exit(yieldSize(bf) <- yieldSize)
        on.exit(close(bf), add=TRUE)
        if (!isOpen(bf))
            bf <- open(bf)
        ## here we do not do 'isDuplicate=FALSE' because it may slow down
        ## the scanning by an order of magnitude
        sbf <- scanBamFlag(isPaired=!singleEnd,
                           isUnmappedQuery=FALSE,
                           isSecondaryAlignment=FALSE,
                           isNotPassingQualityControls=FALSE)
        w <- scanBam(bf, param=ScanBamParam(flag=sbf, what="qwidth"))[[1]]$qwidth
        if (length(w) > 0) {
          tab <- table(w)
          w <- as.integer(names(tab[which.max(tab)]))
        } else {
          w <- 0
          warning(sprintf(paste("Could not figure out read length from",
                                "BAM file %s. Possibly something wrong with it",
                                sep="\n"), basename(bf)))
        }
        w
    }
    if (length(bfl) > 1 && bpnworkers(BPPARAM) > 1) {
        bpparsilent <- BPPARAM
        bpprogressbar(bpparsilent) <- FALSE
        rlen <- unlist(bplapply(bfl, qw, singleEnd=singleEnd,
                                BPPARAM=bpparsilent),
                       use.names=FALSE)
    } else
        rlen <- unname(vapply(bfl, qw, singleEnd=singleEnd,
                              FUN.VALUE=integer(1L)))
    ## sometimes shorter read lengths arise from read trimming, we identify
    ## those shorter lengths as less frequent in the whole dataset and replace
    ## them by the most frequent one, which should be the original layout
    freq <- table(rlen) / length(rlen)
    highfreqrlen <- as.integer(names(freq)[which.max(freq)])
    lowfreqrlen <- as.integer(names(freq)[freq < 0.1])
    rlen[rlen %in% lowfreqrlen] <- highfreqrlen
    rlen <- as.integer(rlen)

    if (verbose) {
      urlen <- sort(unique(rlen))
      rlenstr <- sprintf("%d (%s%dnt)", table(rlen),
                         ifelse(singleEnd, "", "2x"), urlen)
      rlenstr <- sprintf("%s, %s",
                         ifelse(singleEnd, "single-end",
                                "paired-end"),
                         paste(rlenstr, collapse=", "))
      cli_alert_info("Library layout: {rlenstr}.")
    }

    rlen
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
#' @importFrom cli cli_alert_success
.loadAnnotationPackageObject <- function(pkgName, argName, pkgType,
                                         verbose=TRUE) {

    callobj <- match.call()
    annotObj <- NULL

    if (is.character(pkgName)) {
        if (!pkgName %in% installed.packages(noCache=TRUE)[, "Package"])
            stop(sprintf("Please install the Bioconductor package %s.", 
                         pkgName))
        if (!.isPkgLoaded(pkgName)) {
            loaded <- suppressPackageStartupMessages(requireNamespace(pkgName))
            if (loaded && verbose)
                cli_alert_success("Loaded {pkgType} annotation package {pkgName}.")
            else if (!loaded)
                stop(sprintf("Package %s could not be loaded.", pkgName))
        }
        tryCatch({
            if (!paste0("package:", pkgName) %in% search())
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
#' @importFrom GenomeInfoDb seqlevelsStyle<- seqinfo seqinfo<- seqlevels
#' @importFrom GenomeInfoDb genome genome<-
.matchSeqinfo <- function(gal, tx, verbose=TRUE) {
    stopifnot("GAlignments" %in% class(gal) ||
              "GAlignmentPairs" %in% class(gal) ||
              "GAlignmentsList" %in% class(gal) ||
              "TxDb" %in% class(tx)) ## QC

    if (length(intersect(seqlevelsStyle(gal), seqlevelsStyle(tx))) > 0)
      return(gal)

    if (is.na(genome(gal)[1]))
      genome(gal) <- genome(tx)

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
            stop("None of the chromosomes in the input BAM file have the ",
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

## private function .checkSEandSCargs()
.checkSEandSCargs <- function(singleEnd, stdChrom) {
    if (!missing(singleEnd)) {
        if (!is.logical(singleEnd) || length(singleEnd) != 1)
            stop("'singleEnd' must be a single logical value")
    }
    if (!missing(stdChrom)) {
        if (!is.logical(stdChrom) || length(stdChrom) != 1)
            stop("'stdChrom' must be a logical value")
    }
}


## private function .estimateStrandedness() it assumes that .checkBamFiles(),
## .checkYieldSize(), .checkPairedEnd() and .checkBamFileListArgs() have
## been called before
#' @importFrom methods is
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_progress_bar cli_progress_done
.estimateStrandedness <- function(bfl, txdb, singleEnd, stdChrom, exonsBy,
                                  minnaln, verbose,
                                  BPPARAM=SerialParam(progressbar=verbose)) {
    stopifnot(is(bfl, "BamFileList")) ## QC

    if (is.character(txdb))
        txdb <- .loadAnnotationPackageObject(txdb, "txdb", "TxDb",
                                             verbose=verbose)

    annot <- exonsBy(txdb, by=exonsBy)
    if (stdChrom) {
        annot <- keepStandardChromosomes(annot, pruning.mode="fine")
        annot <- annot[lengths(annot) > 0]
    }
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE,
                           isProperPair=!singleEnd,
                           isSecondaryAlignment=FALSE,
                           isNotPassingQualityControls=FALSE)
    param <- ScanBamParam(flag=sbflags)

    if (length(bfl) > 1 && bpnworkers(BPPARAM) > 1) {
        verbose <- FALSE
        strbysm <- bplapply(bfl, .strness_oneBAM, tx=annot, stdChrom=stdChrom,
                            singleEnd=singleEnd, strandMode=1L, param=param,
                            minnaln=minnaln, verbose=verbose, BPPARAM=BPPARAM)
    } else {
        if (verbose)
            idpb <- cli_progress_bar("Estimating strandedness", total=length(bfl))
        strbysm <- lapply(bfl, .strness_oneBAM, tx=annot, stdChrom=stdChrom,
                          singleEnd=singleEnd, strandMode=1L, param=param,
                          minnaln=minnaln, verbose=verbose, idpb=idpb)
        if (verbose)
            cli_progress_done(idpb)
    }
    names(strbysm) <- gsub(pattern = ".bam", "", names(strbysm), fixed = TRUE)
    strbysm <- do.call("rbind", strbysm)
    .checkMinNaln(strbysm, minnaln)
    as.data.frame(strbysm)
}
