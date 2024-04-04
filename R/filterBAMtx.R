#' Filter alignments in a BAM file using a transcriptome
#'
#' Filter alignments in a BAM file using criteria based on a
#' transcriptome annotation.
#'
#' @param object gDNAx object obtained with the function 'gDNAdx()'.
#'
#' @param path Directory where to write the output BAM files.
#'
#' @param txflag A value from a call to the function 'filterBAMtxFlag()'.
#'
#' @param param A 'ScanBamParam' object.
#'
#' @param yieldSize (Default 1e6) Number of records in the input BAM file to
#' yield each time the file is read. The lower the value, the smaller memory
#' consumption, but in the case of large BAM files, values below 1e6 records
#' may decrease the overall performance.
#'
#' @param wsize (Default 1000) Window size employed when the argument
#' \code{txflag} includes the value \code{isInStrandedWindow=TRUE}.
#'
#' @param wstep (Default 100) Window step employed when the argument
#' \code{txflag} includes the value \code{isInStrandedWindow=TRUE}.
#'
#' @param pstrness (Default 0.6) Strandedness value above which we consider
#' a target read alignment to occur in stranded window.
#'
#' @param p.value (Default 0.05) Numeric value between 0 and 1 specifying the
#' adjusted p-value cutoff under which we reject the null hypothesis that a
#' target read alignment occurs in a window with an strandedness value below
#' the one given in the parameter \code{pstrness}. This parameter is only used
#' when the argument \code{txflag} includes the value
#' \code{isInStrandedWindow=TRUE}.
#'
#' @param p.adj.method (Default "holm") Method used to adjust p-values that
#' are compared against the cutoff value specified in the parameter
#' \code{p.value}. Adjusted p-values are calculated using the base R function
#' \code{p.adjust()} and this argument is directly passed to the argument
#' \code{method} of that function.
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
#' @return A vector of output filename paths.
#' 
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Getting the 'gDNAx' object
#' bamfiles <- LiYu22subsetBAMfiles()
#' bamfiles <- bamfiles[c(1,7)] # using a subset of samples
#' gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' 
#' # Filtering splice-compatible alignments and writing them into new BAM files
#' fbf <- filterBAMtxFlag(isSpliceCompatibleJunction=TRUE,
#'                        isSpliceCompatibleExonic=TRUE)
#' dir <- tempdir()
#' fstats <- filterBAMtx(gdnax, path=dir, txflag=fbf)
#' list.files(dir, pattern="*.bam$")
#' 
#'
#' @importFrom S4Vectors mcols
#' @importFrom Rsamtools BamFileList scanBamFlag ScanBamParam
#' @importFrom Rsamtools bamWhat bamWhat<-
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom BiocParallel SerialParam bplapply bpnworkers multicoreWorkers
#' @importFrom cli cli_alert_danger cli_alert_warning cli_progress_bar cli_progress_done
#' @importFrom stats p.adjust.methods
#' @export
#' @rdname filterBAMtx
filterBAMtx <- function(object, path=".", txflag=filterBAMtxFlag(),
                        param=ScanBamParam(), yieldSize=1000000,
                        wsize=1000, wstep=100, pstrness=0.6, p.value=0.05,
                        p.adj.method=p.adjust.methods,
                        verbose=TRUE, BPPARAM=SerialParam(progressbar=verbose)) {
    bfl <- object@bfl
    singleEnd <- object@singleEnd
    strandMode <- object@strandMode
    stdChrom <- object@stdChrom
    igc <- object@intergenic
    int <- object@intronic
    tx <- object@transcripts
    tx2gene <- object@tx2gene

    if (!file.exists(path))
        stop(sprintf("path %s does not exist.", path))
    p.adj.method <- match.arg(p.adj.method)

    if (txflag == 0)
        stop("No alignment type selected in argument 'txflag'. Please use ",
             "the function 'filterBAMtxFlag()' to select at least one.")
    if (testBAMtxFlag(txflag, "isInStrandedWindow") && is.na(strandMode) &&
        verbose) {
        cli_alert_danger("Data is unstranded, stranded windows will not help.")
    }

    yieldSize <- .checkYieldSize(yieldSize)
    bfl <- lapply(bfl, function(x, ys) {
                            yieldSize(x) <- ys
                            x
                        }, yieldSize)
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    what0 <- c("rname", "strand", "pos", "cigar", "qname")
    if (singleEnd)
        bamWhat(param) <- setdiff(bamWhat(param), c("groupid", "mate_status"))
    else {
        flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE,
                             isUnmappedQuery=FALSE)
        what0 <- c(what0, "flag", "groupid", "mate_status")
    }
    param <- GenomicAlignments:::.normargParam(param, flag0, what0)

    if (verbose)
        cli_alert_info(sprintf("gDNAx object for %d BAM files", length(bfl)))
    if (verbose && length(bfl) > 10 && bpnworkers(BPPARAM) == 1 &&
        multicoreWorkers() > 1) {
        fmtstr <- "%d BAM files, consider setting the 'BPPARAM' argument."
        cli_alert_warning(sprintf(fmtstr, length(bfl)))
    }

    out.st <- NULL
    if (length(bfl) > 1 && bpnworkers(BPPARAM) > 1) {
        verbose <- FALSE
        out.st <- bplapply(bfl, .filter_oneBAMtx, igc=igc, int=int, tx=tx,
                           path=path, txflag=txflag, singleEnd=singleEnd,
                           strandMode=strandMode, stdChrom=stdChrom,
                           tx2gene=tx2gene, param=param, wsize=wsize,
                           wstep=wstep, pstrness=pstrness, p.value=p.value,
                           p.adj.method=p.adjust.methods,
                           verbose=verbose, BPPARAM=BPPARAM)
    } else {
      idpb <- NULL
      if (verbose)
          idpb <- cli_progress_bar("Filtering", total=length(bfl))
      out.st <- lapply(bfl, .filter_oneBAMtx, igc=igc, int=int,tx=tx,path=path,
                       txflag=txflag, singleEnd=singleEnd,
                       strandMode=strandMode, stdChrom=stdChrom,
                       tx2gene=tx2gene, param=param, wsize=wsize, wstep=wstep,
                       pstrness=pstrness, p.value=p.value,
                       p.adj.method=p.adjust.methods, verbose=verbose,
                       idpb=idpb)
      if (verbose)
          cli_progress_done(idpb)
    }

    out.st <- data.frame(do.call("rbind", out.st))
    rownames(out.st) <- gsub(".bam", "", rownames(out.st))

    out.st
}

#' @importFrom BiocGenerics basename path
#' @importFrom S4Vectors FilterRules
#' @importFrom Rsamtools filterBam sortBam indexBam
#' @importFrom cli cli_progress_update
.filter_oneBAMtx <- function(bf, igc, int, tx, path, txflag, singleEnd,
                             strandMode, stdChrom, tx2gene, param, wsize, wstep,
                             pstrness, p.value, p.adj.method, verbose, idpb) {

    onesuffix <- c(isIntergenic="IGC",
                   isIntronic="INT",
                   isSpliceCompatibleJunction="SCJ",
                   isSpliceCompatibleExonic="SCE",
                   isInStrandedWindow="STW")
    suffix <- "_"
    for (flag in TXFLAG_BITNAMES)
        if (testBAMtxFlag(txflag, flag))
            suffix <- paste0(suffix, onesuffix[flag])

    bamoutfile <- sprintf("%s/%s.bam", path,
                        paste0(gsub(".bam", "", basename(path(bf))), suffix))
    baioutfile <- sprintf("%s/%s.bai", path,
                        paste0(gsub(".bam", "", basename(path(bf))), suffix))
    statsenvname <- sprintf("stats_%s", gsub(".bam", "", basename(path(bf))))
    assign(statsenvname, new.env())
    assign("stats", c(NALN=0L, NIGC=0L, NINT=0L, NSCJ=0L, NSCE=0L, NSTW=0L,
                      NNCH=0L), envir=get(statsenvname))

    filter <- FilterRules(list(BAMtx=.bamtx_filter))
    ff <- filterBam(bf, bamoutfile, param=param, filter=filter,
                    indexDestination=FALSE)
    if (!singleEnd) {
      ## need to do this on the same directory of the output BAM file
      ## to avoid breaking Rsamtools::filterBam() with 'indexDestination=TRUE'
      ## when /tmp is in a different partition as the output BAM file
      ## leading to the error "cannot rename file .. invalid cross-device link"
      tmpfname <- basename(tempfile())
      sortBam(bamoutfile, file.path(dirname(bamoutfile), tmpfname))
      file.rename(file.path(dirname(bamoutfile), paste0(tmpfname, ".bam")),
                            bamoutfile)
    }
    baifname <- indexBam(bamoutfile)
    file.rename(baifname, baioutfile)
    stats <- get("stats", envir=get(statsenvname))
    for (flag in TXFLAG_BITNAMES)
        if (!testBAMtxFlag(txflag, flag))
            stats[paste0("N", onesuffix[flag])] <- NA_integer_

    if (verbose)
        cli_progress_update(id=idpb)

    stats
}

#' @importFrom Rsamtools bamWhat bamTag
#' @importFrom S4Vectors DataFrame mcols mcols<-
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicAlignments GAlignments njunc first
.bamtx_filter <- function(x) {

    n <- 5  ## this number is derived from the fact that .bamtx_filter()
            ## is called by 'eval()' within the 'filterBam()' function
            ## and allows one to access the objects in the scope of
            ## .filter_oneBAMtx() through the environment 'parent.frame(n)'
    bf <- get("bf", envir=parent.frame(n))
    param <- get("param", envir=parent.frame(n))
    txflag <- get("txflag", envir=parent.frame(n))
    singleEnd <- get("singleEnd", envir=parent.frame(n))
    strandMode <- get("strandMode", envir=parent.frame(n))
    stdChrom <- get("stdChrom", envir=parent.frame(n))
    igc <- get("igc", envir=parent.frame(n))
    int <- get("int", envir=parent.frame(n))
    tx <- get("tx", envir=parent.frame(n))
    tx2gene <- get("tx2gene", envir=parent.frame(n))
    verbose <- get("verbose", envir=parent.frame(n))
    wsize <- get("wsize", envir=parent.frame(n))
    wstep <- get("wstep", envir=parent.frame(n))
    pstrness <- get("pstrness", envir=parent.frame(n))
    p.value <- get("p.value", envir=parent.frame(n))
    p.adj.method <- get("p.adj.method", envir=parent.frame(n))
    statsenvname <- get("statsenvname", envir=parent.frame(n))
    statsenv <- get(statsenvname, envir=parent.frame(n))

    seqlengths <- .getSeqlen(bf, x)
    gal <- GAlignments(seqnames=x$rname, pos=x$pos, cigar=x$cigar,
                       strand=x$strand, seqlengths=seqlengths)
    stopifnot(nrow(x) == length(gal)) ## QC
    cnames <- setdiff(c(bamWhat(param), bamTag(param)),
                        c("rname", "pos", "cigar", "strand"))
    dtf <- DataFrame(idx=1:nrow(x))
    if (length(cnames) > 0) {
        dtf2 <- do.call(DataFrame, as.list(x[cnames]))
        colnames(dtf2) <- cnames
        dtf <- cbind(dtf2, dtf)
    }
    mcols(gal) <- dtf
    if (!singleEnd)
        gal <- .makeGALPE(param, strandMode, gal)
    envstats <- get("stats", envir=statsenv)
    if (stdChrom) {
        ngal <- length(gal)
        if (singleEnd)
            gal <- keepStandardChromosomes(gal, pruning.mode="coarse")
        else
            gal <- keepStandardChromosomes(gal, pruning.mode="fine")
        if (length(gal) < ngal)
          envstats["NNCH"] <- envstats["NNCH"] + ngal - length(gal)
        if (length(gal) == 0) { ## if all alignments are in non-standard chromosomes
          assign("stats", envstats, envir=statsenv)
          return(rep(FALSE, nrow(x)))
        }
    }
    
    gal <- .matchSeqinfo(gal, tx, verbose)
    
    alntype <- .getalntype(gal, txflag, igc, int, strandMode, tx, tx2gene,
                           singleEnd, wsize, wstep, pstrness, p.value,
                           p.adj.method)
    envstats <- envstats + alntype$stats
    assign("stats", envstats, envir=statsenv)
    
    mask <- rep(FALSE, nrow(x))
    if (any(alntype$mask)) {
        keep <- integer(0)
        if (singleEnd)
            keep <- mcols(gal)$idx[alntype$mask]
        else
            keep <- c(mcols(first(gal))$idx[alntype$mask],
                      mcols(last(gal))$idx[alntype$mask])
        mask[keep] <- TRUE
    }
    mask
}

## private function .getSeqlen()
#' @importFrom BiocGenerics path
#' @importFrom GenomeInfoDb seqlengths
.getSeqlen <- function(bf, x) {
    seqlengths <- seqlengths(bf)
    if (!is.null(seqlengths)) {
        bad <- setdiff(levels(x$rname), names(seqlengths))
        if (length(bad) > 0) {
            bad <- paste(bad, collapse="' '")
            msg <- sprintf(paste("'rname' lengths not in BamFile header;",
                                 "seqlengths not used\n  file: %s\n  missing",
                                 "rname(s): '%s'", sep=" "), path(bf), bad)
          warning(msg)
          seqlengths <- NULL
        }
    }
    return(seqlengths)
}

## private function .makeGALPE()
#' @importFrom Rsamtools bamWhat bamTag
#' @importFrom S4Vectors mcols
.makeGALPE <- function(param, strandMode, gal) {
    use.mcols <- setdiff(c(bamWhat(param), bamTag(param), colnames(mcols(gal))),
                         c("rname", "pos", "cigar", "strand", "qname", "flag",
                           "groupid", "mate_status"))
    strandMode2 <- strandMode
    if (is.na(strandMode))
        strandMode2 <- 1L
    makeGALP <- GenomicAlignments:::.make_GAlignmentPairs_from_GAlignments
    suppressWarnings(gal <- makeGALP(gal, strandMode2, use.mcols=use.mcols))
    return(gal)
}

## private function .getalntype()
.getalntype <- function(gal, txflag, igc, int, strandMode, tx, tx2gene,
                        singleEnd, wsize, wstep, pstrness, p.value,
                        p.adj.method) {

    whalnstr <- character(0)
    stats <- c(NALN=length(gal), NIGC=0L, NINT=0L, NSCJ=0L, NSCE=0L, NSTW=0L,
               NNCH=0L)
    if (length(gal) == 0)
      return(list(mask=logical(0), whalnstr=whalnstr, stats=stats))

    tmask <- imask <- rep(FALSE, length(gal))
    if (testBAMtxFlag(txflag, "isIntergenic")) {
        igcaln <- .igcAlignments(gal, igc, fragmentsLen=FALSE)
        tmask <- tmask | igcaln$igcmask
        whalnstr <- c(whalnstr, "IGC")
        stats["NIGC"] <- sum(igcaln$igcmask)
    }
    if (testBAMtxFlag(txflag, "isIntronic")) {
        intaln <- .intAlignments(gal, int, strandMode, fragmentsLen=FALSE)
        tmask <- tmask | intaln$intmask
        whalnstr <- c(whalnstr, "INT")
        stats["NINT"] <- sum(intaln$intmask)
    }
    if (testBAMtxFlag(txflag, "isSpliceCompatibleJunction") ||
        testBAMtxFlag(txflag, "isSpliceCompatibleExonic")) {
        scoaln <- .scoAlignments(gal, tx, tx2gene, singleEnd, strandMode,
                                 fragmentsLen=FALSE)
        if (testBAMtxFlag(txflag, "isSpliceCompatibleJunction")) {
            imask <- scoaln$scjmask
            tmask <- tmask | scoaln$scjmask
            whalnstr <- c(whalnstr, "SCJ")
            stats["NSCJ"] <- sum(scoaln$scjmask)
        }
        if (testBAMtxFlag(txflag, "isSpliceCompatibleExonic")) {
            tmask <- tmask | scoaln$scemask
            whalnstr <- c(whalnstr, "SCE")
            stats["NSCE"] <- sum(scoaln$scemask)
        }
    }
    if (testBAMtxFlag(txflag, "isInStrandedWindow")) {
        stwmask <- .strandedWindowMask(gal, tx, singleEnd, wsize, wstep,
                                       pstrness, p.value, p.adj.method, tmask,
                                       imask)
        tmask <- tmask & stwmask
        whalnstr <- c(whalnstr, "STW")
        stats["NSTW"] <- sum(stwmask) 
    }
    return(list(mask=tmask, whalnstr=whalnstr, stats=stats))
}


TXFLAG_BITNAMES <- c("isIntergenic",
                     "isIntronic",
                     "isSpliceCompatibleJunction",
                     "isSpliceCompatibleExonic",
                     "isInStrandedWindow")

## adapted from Rsamtools::scanBamFlag()

#' Transcriptome-based parameters for filtering BAM files
#'
#' Use 'filterBAMtxFlag()' to set what types of alignment in a BAM 
#' file should be filtered using the function 'filterBAMtx()',
#' among being splice-compatible with one or more exon-exon junctions,
#' splice-compatible exonic, splice-compatible exonic in a window, intronic or
#' intergenic.
#' 
#' @param isSpliceCompatibleJunction (Default FALSE) Logical value indicating
#'        if spliced alignments overlapping a transcript in a 
#'        "splice compatible" way should be included in the BAM file. For
#'        paired-end reads, one or both alignments must have one or more splice
#'        sites compatible with splicing. See 
#'        \code{\link[GenomicAlignments:OverlapEncodings-class]{OverlapEncodings}}.
#' 
#' @param isSpliceCompatibleExonic (Default FALSE) Logical value indicating
#'        if alignments without a splice site, but that overlap a transcript
#'        in a "splice compatible" way, should be included in the BAM file.
#'        For paired-end reads, none of the alignments must be spliced, and
#'        each pair can be in different exons (or in the same one), as long as
#'        they are "splice compatible". See 
#'        \code{\link[GenomicAlignments:OverlapEncodings-class]{OverlapEncodings}}.
#'        
#' @param isInStrandedWindow (Default FALSE) Logical value
#'        indicating whether an alignment occurs in a stranded windows. More
#'        concretely, for each alignment, strandedness will be tested with
#'        respect to the rest of the alignments occurring in overlapping
#'        windows. This filter will assign a TRUE value when those tests are
#'        passed.
#'
#' @param isIntronic (Default FALSE) Logical value indicating if alignments
#'        mapping to introns should be included in the BAM file.
#'
#' @param isIntergenic (Default FALSE) Logical value indicating if alignments
#'        aligned to intergenic regions should be included in the BAM file.
#'
#'
#' @examples
#' # Filtering splice-compatible alignments and writing them into new BAM files
#' fbf <- filterBAMtxFlag(isSpliceCompatibleJunction=FALSE,
#'                        isSpliceCompatibleExonic=FALSE,
#'                        isInStrandedWindow=FALSE,
#'                        isIntronic=FALSE,
#'                        isIntergenic=FALSE)
#' 
#' @export
#' @rdname filterBAMtx
filterBAMtxFlag <- function(isSpliceCompatibleJunction=FALSE,
                            isSpliceCompatibleExonic=FALSE,
                            isInStrandedWindow=FALSE,
                            isIntronic=FALSE,
                            isIntergenic=FALSE) {
    flag <- S4Vectors:::makePowersOfTwo(length(TXFLAG_BITNAMES))
    names(flag) <- TXFLAG_BITNAMES
    args <- lapply(as.list(match.call())[-1], eval, parent.frame())
    if (any(vapply(args, length, FUN.VALUE = integer(1L)) > 1L))               
        stop("all arguments must be logical(1)")

    if (length(args) == 0)
        args <- formals(filterBAMtxFlag)

    idx <- names(args[vapply(args, function(x) !is.na(x) && x,
                             FUN.VALUE = logical(1L))])
    keep <- Reduce("+", flag[names(flag) %in% idx], 0L)

    keep
}


#' @param flag A value from a call to the function 'filterBAMtxFlag()'.
#' 
#' @param value A character vector with the name of a flag.
#' 
#' @importFrom bitops bitAnd
#'
#' @examples
#' testBAMtxFlag(fbf, "isSpliceCompatibleJunction")
#'                        
#' @export
#' @rdname filterBAMtx
testBAMtxFlag <- function(flag, value) {
    if (length(value) != 1 || !value %in% TXFLAG_BITNAMES) {
        msg <- sprintf("'is' must be character(1) in '%s'",
                        paste(TXFLAG_BITNAMES, collapse="' '"))
        stop(msg)
    }
    i <- 2^(match(value, TXFLAG_BITNAMES) - 1L)
    bitAnd(flag, i) == i
}

## private function .strandedWindowMask() to identify alignments that are
## stranded using a windowing approach
#' @importFrom stats p.adjust pbinom
#' @importFrom S4Vectors split decode pc
#' @importFrom IRanges coverage resize width slidingWindows Views viewSums
#' @importFrom IRanges RleList IntegerList relist lapply
#' @importFrom GenomicRanges GRanges GRangesList grglist granges ranges trim
#' @importFrom GenomicRanges seqnames start start<- end end<- match
#' @importFrom GenomicAlignments first last
#' @importFrom matrixStats rowMins

.strandedWindowMask <- function(gal, tx, singleEnd, wsize, wstep, pstrness,
                                p.value, p.adj.method, tmask, imask) {

    if (length(gal) == 0)
      return(logical(0))

    ## generate coverage separately for the forward and reverse strands
    grlfwd <- grlrev <- GRangesList()
    if (!singleEnd) {
        ## separately for each mate in the case of paired-end read alignments
        gal1 <- first(gal, real.strand=TRUE)
        gal2 <- last(gal, real.strand=TRUE)
        grl1fwd <- grglist(gal1[strand(gal1) == "+"])
        grl1rev <- grglist(gal1[strand(gal1) == "-"])
        grl2fwd <- grglist(gal2[strand(gal2) == "+"])
        grl2rev <- grglist(gal2[strand(gal2) == "-"])
        grlfwd <- sort(c(unlist(grl1fwd), unlist(grl2fwd)))
        grlrev <- sort(c(unlist(grl1rev), unlist(grl2rev)))
    } else {
        grlfwd <- grglist(gal[strand(gal) == "+"])
        grlrev <- grglist(gal[strand(gal) == "-"])
    }

    covfwd <- coverage(grlfwd)
    covrev <- coverage(grlrev)

    ## fetch genomic ranges for reads where we are going to test for
    ## strandedness, restricting calculations to those that may have
    ## been already previously selected through the input parameter 'tmask'.
    ## if no alignment is selected in 'tmask', then do calculations in all
    ## of the alignments. alignments 'imask' are selected to be included
    ## anyway, these can be also excluded from the calculations
    galtarget <- gal
    tmask <- tmask & !imask
    if (sum(tmask) > 0)
      galtarget <- gal[tmask]

    gr <- GRanges()
    if (!singleEnd) {
      gr1 <- granges(first(galtarget, real.strand=TRUE))
      gr2 <- granges(last(galtarget, real.strand=TRUE))
      stopifnot(identical(strand(gr1), strand(gr2))) ## QC
      gr <- sort(c(gr1, gr2))
    } else
      gr <- granges(galtarget)

    ## fetch unique genomic ranges
    uniqgr <- unique(gr)

    ## create windows anchored at the ends of each read alignment in 'galtarget'
    suppressWarnings(uniqgrwin <- resize(uniqgr, width(uniqgr)+wsize*2,
                                         fix="center"))
    ## set proper bounds if necessary
    whoutofbounds <- which(trim(uniqgrwin) != uniqgrwin)
    if (length(whoutofbounds) > 0) {
        ## this will result in a smaller initial window range, leading to fewer
        ## sliding windows per read alignment target, i.e., not all read alignment
        ## targets will have the same number of sliding windows ('nwin' below)
        
        mask <- start(uniqgrwin) < 1
        suppressWarnings(start(uniqgrwin[mask]) <- 1)
        sl <- seqlengths(uniqgrwin)[decode(seqnames(uniqgrwin))]
        mask <- end(uniqgrwin) > sl
        suppressWarnings(end(uniqgrwin[mask]) <- sl[mask])
    }

    slw <- slidingWindows(uniqgrwin, width=wsize, step=wstep)
    nwin <- unique(lengths(slw))
    if (length(nwin) > 1) {
        ## if there are different number of windows per read alignment target,
        ## pad those collections of windows with 0-width windows to enforce all
        ## having the maximum possible number of windows
        nwin <- lengths(slw)
        maxnwin <- max(nwin)
        slw2 <- slw[nwin < maxnwin]
        n2pad <- maxnwin - lengths(slw2)
        seqs2pad <- unlist(seqnames(slw2), use.names=FALSE)[cumsum(lengths(slw2))]
        strands2pad <- unlist(strand(slw2), use.names=FALSE)[cumsum(lengths(slw2))]
        gr2pad <- GRanges(rep(seqs2pad, times=n2pad), IRanges(1, 0),
                          strand=rep(strands2pad, times=n2pad))
        gr2pad <- split(gr2pad, rep(1:length(slw2), times=n2pad))
        slw2 <- pc(gr2pad, slw2)
        slw[nwin < maxnwin] <- slw2
        rm(slw2)
        nwin <- unique(lengths(slw))
    }
    stopifnot(length(nwin) == 1) ## QC

    ## sum forward and reverse coverage over windows
    viewsfwd <- Views(covfwd, unlist(slw))
    cntfwd <- viewSums(viewsfwd)
    viewsrev <- Views(covrev, unlist(slw))
    cntrev <- viewSums(viewsrev)

    ## scale counts to library size to avoid counting reads more than once
    scntfwd <- cntfwd
    scntrev <- cntrev
    zeromask <- all(cntfwd == 0L)
    scntfwd[!zeromask] <- ceiling(length(gal) *
                                  (cntfwd[!zeromask] / sum(covfwd[!zeromask])))
    zeromask <- all(cntrev == 0L)
    scntrev[!zeromask] <- ceiling(length(gal) *
                                  (cntrev[!zeromask] / sum(covrev[!zeromask])))
    scnttot <- scntfwd + scntrev

    ## set alignments counts according to their strand: if strand is '+' take
    ## the counts on the '+' strand, if strand is '-', take the counts on the
    ## '-' strand, and if strand is '*' take the maximum value between the
    ## counts in the positive and negative strands
    strbyseq <- split(rep(strand(uniqgr), each=nwin),
                      rep(seqnames(uniqgr), each=nwin))
    scnt <- scntfwd
    strmask <- strbyseq == "-"
    ## if (any(strmask))
    ##   scnt[strmask] <- scntrev[strmask]
    scnt[strmask] <- scntrev[strmask]
    strmask <- strbyseq == "*"
    ## if (any(strmask))
    ##   scnt[strmask] <- pmax(scntfwd[strmask], scntrev[strmask])
    scnt[strmask] <- IntegerList(mapply(pmax, scntfwd[strmask],
                                        scntrev[strmask]))

    ## for each window, calculate the p-value of a one-sided exact binomial test
    ## with a given hypothesized strandedness proportion above which we consider
    ## a target read alignment to occur in a stranded window
    p <- pbinom(unlist(scnt, use.names=FALSE), unlist(scnttot, use.names=FALSE),
                prob=pstrness, lower.tail=FALSE)
    p <- relist(p, scnt)

    ## re-arrange window p-values by unique target read alignment, select the
    ## smallest p-value among the windows of each unique target read alignment,
    ## and adjust them using the Bonferroni method
    p <- lapply(p, matrix, ncol=nwin, byrow=TRUE)
    p <- lapply(p, rowMins)
    p <- lapply(p, p.adjust, method=p.adj.method)

    uniqgr$p <- unlist(p, use.names=FALSE)

    ## match back p-values to the target read alignments. in the case of
    ## paired-end read alignments, select the smallest p-value between the two
    ## mates
    if (!singleEnd) {
      gr1 <- granges(first(galtarget, real.strand=TRUE))
      gr2 <- granges(last(galtarget, real.strand=TRUE))
      mt1 <- match(gr1, uniqgr)
      stopifnot(all(!is.na(mt1))) ## QC
      mt2 <- match(gr2, uniqgr)
      stopifnot(all(!is.na(mt2))) ## QC
      p1 <- uniqgr$p[mt1]
      p2 <- uniqgr$p[mt2]
      p <- pmin(p1, p2)
    } else {
      gr <- granges(galtarget)
      mt <- match(gr, uniqgr)
      p <- uniqgr$p[mt]
    }

    ## build the final logical mask that selects target read alignments with a
    ## p-value below a given threshold
    mask <- p < p.value
    if (sum(tmask) > 0) {
      mask <- rep(FALSE, length(gal))
      mask[tmask] <- p < p.value
    }
    ## include alignments selected in 'imask'
    mask <- mask | imask

    return(mask)
}
