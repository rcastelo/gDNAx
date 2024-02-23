#' Calculate gDNA diagnostics
#'
#' Calculate diagnostics for assessing the presence of genomic DNA (gDNA)
#' in RNA-seq data over a subset of the alignments in the input BAM files.
#'
#' @param bfl A \code{BamFile} or \code{BamFileList} object, or a character
#' string vector of BAM filenames.
#'
#' @param txdb A character string of a \code{TxDb} package, or a \code{TxDb}
#' object, with gene and transcript annotations. For accurate calculations, it
#' is important that the version of these annotations matches the version of
#' the annotations used to inform the alignment of spliced reads, by the
#' short-read aligner software that generated the input BAM files.
#'
#' @param singleEnd (Default FALSE) Logical value indicating if reads are
#' single (\code{TRUE}) or paired-end (\code{FALSE}).
#'
#' @param strandMode (Default 1L) Numeric vector which can take values 0, 1,
#' 2 or \code{NA}. The strand mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail. If \code{singleEnd = TRUE}, then \code{strandMode}
#' is ignored. For not strand-specific libraries, use \code{NA}
#'
#' @param stdChrom (Default TRUE) Logical value indicating whether only
#' alignments in the 'standard chromosomes' should be used. Consult the help
#' page of the function \code{\link[GenomeInfoDb]{keepStandardChromosomes}}
#' from the package \code{GenomeInfoDb} for further
#' information.
#'
#' @param yieldSize (Default 1e5) Number of records to read from each input BAM
#' file to calculate the diagnostics.
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
#' @return A \linkS4class{gDNAx} object.
#' 
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Retrieving BAM files
#' bamfiles <- LiYu22subsetBAMfiles()
#' bamfiles <- bamfiles[c(1,4,7)] # using a subset of samples
#' 
#' # Getting information about the gDNA concentrations of each BAM file
#' pdat <- LiYu22phenoData(bamfiles)
#' 
#' gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' gdnax
#'
#' @importFrom BiocGenerics basename path
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomeInfoDb keepStandardChromosomes genome
#' @importFrom BiocParallel SerialParam bplapply bpnworkers bpprogressbar<-
#' @importFrom methods new
#' @export
#' @rdname gDNAdx
gDNAdx <- function(bfl, txdb, singleEnd=TRUE, strandMode=1L, stdChrom=TRUE,
                   yieldSize=100000L, verbose=TRUE,
                   BPPARAM=SerialParam(progressbar=verbose)) {
    yieldSize <- .checkYieldSize(yieldSize)
    bfl <- .checkBamFileListArgs(bfl, singleEnd, fragments=FALSE, yieldSize)
    strandMode <- .checkStrandMode(strandMode)
    .checkOtherArgs(singleEnd, stdChrom)
    if (is.character(txdb))
        txdb <- .loadAnnotationPackageObject(txdb, "txdb", "TxDb",
                                             verbose=verbose)
    bpparsilent <- BPPARAM
    bpprogressbar(bpparsilent) <- FALSE
    peflag <- .allPairedEnd(bfl, singleEnd, bpparsilent)
    rlens <- .readLengths(bfl, singleEnd)
    maxfrglen <- max(rlens)
    if (peflag)
      maxfrglen <- 2 * peflag
    ## minfrglen <- .minFrgLen(bfl, singleEnd)
    if (verbose)
        message(sprintf("Fetching annotations for %s", genome(txdb)[1]))
    igcintrng <- .fetchIGCandINTrng(txdb, maxfrglen, stdChrom, strandMode)
    exbytx <- exonsBy(txdb, by="tx") ## fetch transcript annotations
    if (stdChrom) {
        exbytx <- keepStandardChromosomes(exbytx, pruning.mode="fine")
        exbytx <- exbytx[lengths(exbytx) > 0]
    }
    genetxid <- suppressMessages(select(txdb, keys=names(exbytx),
                                        columns="GENEID", keytype="TXID"))
    stopifnot(identical(as.integer(genetxid$TXID), seq_len(nrow(genetxid))))
    tx2gene <- rep(NA_character_, nrow(genetxid))
    mask <- !is.na(genetxid$GENEID)
    tx2gene[mask] <- genetxid$GENEID[mask]
    rm(genetxid)
    sbflags <- scanBamFlag(isUnmappedQuery=FALSE, isProperPair=!singleEnd,
                           isSecondaryAlignment=FALSE, isDuplicate=FALSE,
                           isNotPassingQualityControls=FALSE)
    param <- ScanBamParam(flag=sbflags)
    
    if (verbose)
        message("Start processing BAM file(s)")
    dxBAMs <- NULL
    if (length(bfl) > 1 && bpnworkers(BPPARAM) > 1) {
        verbose <- FALSE
        dxBAMs <- bplapply(bfl, .gDNAdx_oneBAM, igc=igcintrng$igcrng,
                           int=igcintrng$intrng, tx=exbytx, tx2gene=tx2gene,
                           stdChrom=stdChrom,singleEnd=singleEnd,
                           strandMode=strandMode, param=param,
                           verbose=verbose, BPPARAM=BPPARAM)
    } else
        dxBAMs <- lapply(bfl, .gDNAdx_oneBAM, igc=igcintrng$igcrng,
                         int=igcintrng$intrng, tx=exbytx, tx2gene=tx2gene,
                         stdChrom=stdChrom, singleEnd=singleEnd,
                         strandMode=strandMode, param=param, verbose=verbose)
    
    dxobj <- .collectDiagn(dxBAMs, bfl, rlens, singleEnd, txdb, strandMode,
                           stdChrom, yieldSize, igcintrng, exbytx, tx2gene,
                           verbose)
    dxobj
}


## private function .gDNAdx_oneBAM() to calculate diagnostics from one BAM file
#' @importFrom methods as
#' @importFrom BiocGenerics basename path
#' @importFrom S4Vectors queryHits countQueryHits countSubjectHits
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom GenomicFiles reduceByYield
#' @importFrom Rsamtools isOpen yieldSize
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs
#' @importFrom GenomicAlignments encodeOverlaps isCompatibleWithSplicing
#' @importFrom GenomicAlignments isCompatibleWithSkippedExons granges
#' @importFrom GenomicRanges GRanges grglist
.gDNAdx_oneBAM <- function(bf, igc, int, tx, tx2gene, stdChrom, singleEnd,
                            strandMode, param, verbose) {
    if (isOpen(bf))
        close(bf)
    
    if (verbose)
        message(sprintf("Reading first %d alignments from %s",
                        yieldSize(bf), basename(path(bf))))
    
    gal <- NULL
    if (singleEnd)
        gal <- readGAlignments(bf, param=param, use.names=FALSE)
    else {
        if (is.na(strandMode))
            gal <- readGAlignmentPairs(bf, param=param, use.names=FALSE)
        else
            gal <- readGAlignmentPairs(bf, param=param, strandMode=strandMode,
                                        use.names=FALSE)
    }
    if (stdChrom)
        if (singleEnd)
            gal <- keepStandardChromosomes(gal, pruning.mode="coarse")
        else
            gal <- keepStandardChromosomes(gal, pruning.mode="fine")

    if (verbose)
        message(sprintf("Processing alignments from %s", basename(path(bf))))

    gal <- .matchSeqinfo(gal, tx, verbose)
    naln <- length(gal)

    ## intergenic alignments
    igcaln <- .igcAlignments(gal, igc,fragmentsLen=TRUE,nfrgs=1000)
    ## intronic alignments
    intaln <- .intAlignments(gal, int, strandMode,fragmentsLen=TRUE,nfrgs=1000)
    ## splice-compatible alignments
    scoaln <- .scoAlignments(gal, tx, tx2gene, singleEnd, strandMode,
                             fragmentsLen=TRUE, nfrgs=1000)
    
    ## Assigning reads simultaneously identified as IGC/INT and SCE/SCJ
    ## as SCE/SCJ
    whigc <- (igcaln$igcmask & scoaln$scjmask) | 
                        (igcaln$igcmask & scoaln$scemask)
    igcaln$igcmask[whigc] <- FALSE

    whint <- (intaln$intmask & scoaln$scjmask) | 
                        (intaln$intmask & scoaln$scemask)
    intaln$intmask[whint] <- FALSE

    dx_oneBAM <- .getdx_oneBAM(igcaln, intaln, scoaln, singleEnd,
                                strandMode, gal, tx, naln)
    dx_oneBAM
}

## private function .getdx_oneBAM() 
.getdx_oneBAM <- function(igcaln, intaln, scoaln, singleEnd, strandMode,
                            gal, tx, naln) {
    nigcaln <- sum(igcaln$igcmask)
    nintaln <- sum(intaln$intmask)
    nscjaln <- sum(scoaln$scjmask)
    nscealn <- sum(scoaln$scemask)
    nsccaln <- sum(scoaln$sccmask)
    
    igcfrglen <- intfrglen <- scjfrglen <- scefrglen <- NA_integer_
    if (!singleEnd) { ## fragments length from paired-end reads
        igcfrglen <- igcaln$igcfrglen
        intfrglen <- intaln$intfrglen
        scjfrglen <- scoaln$scjfrglen
        scefrglen <- scoaln$scefrglen
    }
    
    ## strandedness
    if (is.na(strandMode))
        strness <- NA
    else
        strness <- .getStrandedness(gal, tx, reportAll=FALSE)
    
    return(list(naln=naln, nigcaln=nigcaln, nintaln=nintaln, nscjaln=nscjaln,
                nscealn=nscealn, nsccaln=nsccaln, igcfrglen=igcfrglen,
                intfrglen=intfrglen, scjfrglen=scjfrglen, scefrglen=scefrglen,
                strness=strness))
}

## private function .igcAlignments() to find intergenic alignments
#' @importFrom IRanges findOverlaps
.igcAlignments <- function(gal, igc, fragmentsLen, nfrgs=1000) {
    ovig <- findOverlaps(GRanges(gal), igc, type="within", ignore.strand=TRUE)
    igcmask <- countQueryHits(ovig) > 0
    res <- list(igcmask=igcmask)
    if (fragmentsLen) {
        res$igcfrglen <- width(granges(gal[igcmask]))
        if (length(res$igcfrglen) > nfrgs)
            res$igcfrglen <- sample(res$igcfrglen, nfrgs)
    }

    res
}

## private function .intAlignments() to find intronic alignments
#' @importFrom IRanges findOverlaps
.intAlignments <- function(gal, int, strandMode, fragmentsLen, nfrgs=1000) {
    ignore.strand <- ifelse(is.na(strandMode), TRUE, FALSE)
    ovin <- findOverlaps(GRanges(gal), int, type="within", 
                        ignore.strand=ignore.strand)
    intmask <- countQueryHits(ovin) > 0
    res <- list(intmask=intmask)
    if (fragmentsLen) {
        res$intfrglen <- width(granges(gal[intmask]))
        if (length(res$intfrglen) > nfrgs)
            res$intfrglen <- sample(res$intfrglen, nfrgs)
    }

    res
}

## private function .scoAlignments() to find single-gene splice-compatible
## alignments
#' @importFrom S4Vectors countQueryHits
#' @importFrom Biostrings encoding
#' @importFrom IRanges findOverlaps
#' @importFrom GenomicRanges GRanges grglist
#' @importFrom GenomicAlignments encodeOverlaps isCompatibleWithSplicing
#' @importFrom GenomicAlignments isCompatibleWithSkippedExons
.scoAlignments <- function(gal, tx, tx2gene, singleEnd, strandMode,
                           fragmentsLen, nfrgs=1000) {
    
    ignore.strand <- ifelse(is.na(strandMode), TRUE, FALSE)
    ## calculate overlaps between alignments and transcripts
    ovtx <- findOverlaps(GRanges(gal), tx, ignore.strand=ignore.strand)

    ## build mask to select alignments overlapping a single gene
    ovgenes <- tx2gene[subjectHits(ovtx)]
    gbya <- split(ovgenes, queryHits(ovtx))
    gbya <- sapply(gbya, unique)
    gbya <- lengths(gbya)
    onegenealn <- as.integer(names(gbya)[gbya == 1])
    onegenemaskovtx <- queryHits(ovtx) %in% onegenealn

    ## alignments compatible with splicing over the given transcriptome
    ovtxenc <- encodeOverlaps(grglist(gal), tx, hits=ovtx,
                              flip.query.if.wrong.strand=ignore.strand)
    masksplovtx <- isCompatibleWithSplicing(ovtxenc)
    maskskpovtx <- isCompatibleWithSkippedExons(ovtxenc)
    junencpat <- "2--|--2|3--|--3"
    if (singleEnd)
        junencpat <- "2:|3:"
    maskjunovtx <- grepl(junencpat, encoding(ovtxenc))
    ## logical mask for splice-compatible junction alignments
    scjmaskovtx <- (masksplovtx & maskjunovtx) | (maskskpovtx & maskjunovtx)
    exnencpat <- "1--1"
    if (singleEnd)
        exnencpat <- "1:"
    maskexnovtx <- grepl(exnencpat, encoding(ovtxenc))
    ## logical mask for splice-compatible exonic alignments
    scemaskovtx <- (masksplovtx & maskexnovtx) | (maskskpovtx & maskexnovtx)
    rm(masksplovtx, maskskpovtx, maskjunovtx, maskexnovtx)
    ## splice-compatible where at least one read includes a junction
    ovscjtx <- ovtx[scjmaskovtx & onegenemaskovtx]
    mask <- !duplicated(queryHits(ovscjtx)) ## discard redundant info
    ovscjtx <- ovscjtx[mask]
    scjmask <- countQueryHits(ovscjtx) > 0
    ## splice-compatible where reads align completely within exons
    ovscetx <- ovtx[scemaskovtx & onegenemaskovtx]
    mask <- !duplicated(queryHits(ovscetx)) ## discard redundant info
    ovscetx <- ovscetx[mask]
    scemask <- countQueryHits(ovscetx) > 0
    
    res <- list(scjmask=scjmask, scemask=scemask,
                sccmask=sum(njunc(gal) > 0))
    res <- .getfrglen(singleEnd, fragmentsLen, res, gal, tx, ovscjtx, ovscetx, 
                      strandMode, nfrgs)
    res
}

## private function .getfrglen()
.getfrglen <- function(singleEnd, fragmentsLen, res, gal, tx, ovscjtx, ovscetx, 
                       strandMode, nfrgs) {
    if (!singleEnd && fragmentsLen) { ## for paired-end, get fragments length
        res$scjfrglen <- .sampleFragmentsLength(gal, tx, ovscjtx,
                                                strandMode, nfrgs)
        res$scefrglen <- .sampleFragmentsLength(gal, tx, ovscetx,
                                                strandMode, nfrgs)
    }
    res
}

## private function .sampleFragmentsLength() to calculate fragments length in
## transcript coordinates for a sample of transcript-overlapping alignments.
## The sampling strategy is sufficient to estimate the fragments length
## distribution and necessary to minimize the memory footprint of the
## pmapToTranscripts() function.
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom GenomicRanges GRanges start width strand resize
#' @importFrom GenomicFeatures pmapToTranscripts
#' @importFrom GenomicAlignments granges
.sampleFragmentsLength <- function(gal, tx, alnhits, strandMode, nfrgs) {
    ignore.strand <- ifelse(is.na(strandMode), TRUE, FALSE)
    alnhits.sam <- alnhits
    if (length(alnhits.sam) > nfrgs) {
        sam <- sample(seq_along(alnhits), size=nfrgs)
        alnhits.sam <- alnhits[sort(sam)] ## alnhits is SortedByQueryHits
    }
    grsaln <- granges(gal[queryHits(alnhits.sam)])
    txsaln <- tx[subjectHits(alnhits.sam)]
    grsaln.txstart <- pmapToTranscripts(resize(grsaln, width=1, fix="start"), 
                                        txsaln, ignore.strand=ignore.strand)
    grsaln.txend <- pmapToTranscripts(resize(grsaln, width=1, fix="end"), 
                                        txsaln, ignore.strand=ignore.strand)
    stopifnot(all(width(grsaln.txstart) > 0)) ## QC
    stopifnot(all(width(grsaln.txend) > 0)) ## QC
    stopifnot(all(!duplicated(queryHits(alnhits.sam)))) ## QC
    w <- start(grsaln.txend) - start(grsaln.txstart)
    
    # when 'ignore.strand = TRUE' in pmapToTranscripts(), and strand of read 
    # is "-", 'grsaln.txstart' > 'grsaln.txend', resulting in a 'w' < 0
    if (ignore.strand) {
        wh <- which(as.character(strand(grsaln)) == "-")
        w[wh] <- abs(w[wh])
    }
    stopifnot(all(w > 0)) ## QC
    w
}

## private function .fetchIGCandINTrng()
#' @importFrom IRanges reduce
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomeInfoDb keepStandardChromosomes genome
#' @importFrom GenomicRanges strand strand<- gaps intersect width setdiff
.fetchIGCandINTrng <- function(txdb, maxfrglen, stdChrom, strandMode) {
    ## fetch ranges of annotated genes, discard gene information and project
    ## them into genomic coordinates. use exonsBy() instead of genes() to
    ## catch annotations of genes mapping to multiple sequences and including
    ## exons on different strands.
    exbygn <- exonsBy(txdb, by="gene")
    if (stdChrom) ## important to use 'pruning.mode="fine"' here
        exbygn <- keepStandardChromosomes(exbygn, pruning.mode="fine")
    exbygnrng <- unlist(range(exbygn), use.names=FALSE)
    mcols(exbygnrng) <- NULL
    exbygnrng <- reduce(exbygnrng)

    ## fetch ranges of RepeatMasker annot. and project them into genomic coord.
    rmskdb <- .fetchRmsk(txdb)
    ## strand(rmskdb) <- "*"
    if (stdChrom)
        rmskdb <- keepStandardChromosomes(rmskdb, pruning.mode="coarse")
    rmskrng <- reduce(rmskdb)

    ## fetch ranges of projected intron coordinates, i.e. non-overlapping introns
    # pexbygn <- sort(unlist(reduce(exbygn), use.names=FALSE))
    # pintrons <- sort(intersect(gaps(pexbygn), exbygnrng))
    gnrmskrng <- sort(c(unlist(reduce(exbygn), use.names=FALSE), rmskrng))
    pintrons <- sort(intersect(gaps(gnrmskrng), exbygnrng))
    if (is.na(strandMode))
        pintrons <- setdiff(pintrons, gnrmskrng, ignore.strand=TRUE)
    if (stdChrom)
        pintrons <- keepStandardChromosomes(pintrons, pruning.mode="coarse")
    pintrons <- pintrons[width(pintrons) >= maxfrglen]

    ## put together gene and repeat annotations
    rmskrng <- reduce(rmskdb, ignore.strand=TRUE)
    strand(exbygnrng) <- "*"
    gnrmskrng <- sort(c(exbygnrng, rmskrng))

    ## fetch ranges of intergenic regions as those complementary to the
    ## previous genic and repeat ranges, collapsing neighboring ranges closer
    ## to each other than the maximum sequencing fragment length and/or below
    ## that size, to avoid generating intergenic regions under such a size
    gnrmskrng <- reduce(gnrmskrng, min.gapwidth=maxfrglen)
    gnrmskrng <- sort(gnrmskrng)
    igcrng <- gaps(gnrmskrng)
    igcrng <- igcrng[strand(igcrng) == "*"]
    igcrng <- sort(igcrng)
    igcrng <- igcrng[width(igcrng) >= maxfrglen]
    if (stdChrom)
        igcrng <- keepStandardChromosomes(igcrng, pruning.mode="coarse")
    
    list(igcrng=igcrng, intrng=pintrons)
}

## private function .fetchRmsk()
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom GenomeInfoDb genome
.fetchRmsk <- function(txdb) {
    suppressMessages(ah <- AnnotationHub())
    suppressMessages(ahres <- query(ah, c("UCSCRepeatMasker",genome(txdb)[1])))
    mt <- gregexpr(pattern=" \\([A-Za-z0-9]+\\) ", ahres$title)
    ahresdates <- substr(ahres$title, unlist(mt)+2,
                            unlist(mt)+vapply(mt, attr, "match.length",
                                            FUN.VALUE = integer(1L)) - 3)
    ahresdates <- as.Date(paste0("1", ahresdates), "%d%b%Y")
    ## in case > 1 RM annotations available, pick the most recent one
    suppressMessages(rmskdb <- ahres[[which.max(ahresdates)]])
    names(rmskdb) <- mcols(rmskdb) <- NULL
    return(rmskdb)
}

## private function .collectDiagn()
#' @importFrom BiocGenerics path
#' @importFrom methods new
.collectDiagn <- function(dxBAMs, bfl, readLength, singleEnd, txdb, strandMode,
                          stdChrom, yieldSize, igcintrng, exbytx, tx2gene,
                          verbose) {
    if (verbose)
        message("Collecting diagnostics")
    igcpct <- vapply(dxBAMs, function(x)  ## intergenic %
        100 * x$nigcaln / x$naln, FUN.VALUE = numeric(1L))
    intpct <- vapply(dxBAMs, function(x)  ## intronic %
        100 * x$nintaln / x$naln, FUN.VALUE = numeric(1L))
    scjpct <- vapply(dxBAMs, function(x)  ## splice-compatible junction %
        100 * x$nscjaln / x$naln, FUN.VALUE = numeric(1L))
    scepct <- vapply(dxBAMs, function(x)  ## splice-compatible exonic %
        100 * x$nscealn / x$naln, FUN.VALUE = numeric(1L))
    sccpct <- vapply(dxBAMs, function(x) ##njunc() splice-compatible junction %
        100 * x$nsccaln / x$naln, FUN.VALUE = numeric(1L))
    igcfrglen.mean <- vapply(dxBAMs, function(x) mean(x$igcfrglen),
                                FUN.VALUE = numeric(1L))
    intfrglen.mean <- vapply(dxBAMs, function(x) mean(x$intfrglen),
                                FUN.VALUE = numeric(1L))
    scjfrglen.mean <- vapply(dxBAMs, function(x) mean(x$scjfrglen),
                                FUN.VALUE = numeric(1L))
    scefrglen.mean <- vapply(dxBAMs, function(x) mean(x$scefrglen),
                                FUN.VALUE = numeric(1L))
    strness <- vapply(dxBAMs, function(x) x$strness, FUN.VALUE = numeric(1L))
    snames <- gsub(".bam", "", names(igcpct))
    if (any(duplicated(snames))) {
        stopifnot(identical(names(bfl),  names(igcpct)))
        snames <- path(bfl)
    }
    dx <- data.frame(IGC=igcpct, INT=intpct, SCJ=scjpct, SCE=scepct,SCC=sccpct,
                    IGCFLM=igcfrglen.mean, SCJFLM=scjfrglen.mean,
                    SCEFLM=scefrglen.mean, INTFLM=intfrglen.mean,
                    STRAND=strness, row.names=snames)
    igcfrglen <- lapply(dxBAMs, function(x) x$igcfrglen)
    intfrglen <- lapply(dxBAMs, function(x) x$intfrglen)
    scjfrglen <- lapply(dxBAMs, function(x) x$scjfrglen)
    scefrglen <- lapply(dxBAMs, function(x) x$scefrglen)
    names(igcfrglen) <- gsub(".bam", "", names(igcfrglen))
    names(intfrglen) <- gsub(".bam", "", names(intfrglen))
    names(scjfrglen) <- gsub(".bam", "", names(scjfrglen))
    names(scefrglen) <- gsub(".bam", "", names(scefrglen))
    new("gDNAx", bfl=bfl, txdbpkg=txdb$packageName, singleEnd=singleEnd,
        strandMode=strandMode, stdChrom=stdChrom, readLength=readLength,
        yieldSize=yieldSize, diagnostics=dx, igcfrglen=igcfrglen,
        intfrglen=intfrglen, scjfrglen=scjfrglen, scefrglen=scefrglen,
        intergenic=igcintrng$igcrng, intronic=igcintrng$intrng,
        transcripts=exbytx, tx2gene=tx2gene)
}

#' Plot diagnostics
#'
#' Plot diagnostics calculated with gDNAdx()
#'
#' @param x A 'gDNAx' object.
#'
#' @param group A string character vector or a factor, with as many values
#' as BAM files analyzed in 'x', whose values define groups among those
#' BAM files.
#'
#' @param labelpoints (Default FALSE) A logical indicator that labels points
#' in those plots where each point represents a BAM file. Labels correspond
#' to the index number of the BAM file in 'x'.
#' 
#' @param ... Named arguments to be passed to \code{\link[base:plot]{plot}}.
#'
#' @importFrom plotrix thigmophobe
#' @importFrom graphics grid text par legend
#' @importFrom stats setNames
#' 
#' @examples
#' # plot gDNA diagnostic measures
#' plot(gdnax, group=pdat$gDNA, pch=19)
#' 
#' @export
#' @rdname gDNAdx
setMethod("plot", signature(x="gDNAx"),
function(x, group=1L, labelpoints=FALSE, ...) {
    if (nrow(getDx(x)) == 0)
        stop("no diagnostics available in the input 'gDNAx' object.")
    dx <- getDx(x)
    grp <- .setColorGrouping(group, dx)
    group <- grp$group
    grpcol <- grp$col
    xlab <- "Intergenic alignments (%)"
    par(oma = c(4,1,1,1), mfrow=c(2, 2), mar=c(5, 5, 2, 2), xpd=FALSE)
    ## IGC x SCJ
    plot(dx$IGC, dx$SCJ, panel.first=grid(), col=grpcol[as.integer(group)],
         xlab=xlab, ylab="Splice-comp. junction aln. (%)", las=1, ...)
    if (labelpoints) {
        pos <- setNames(thigmophobe(dx$IGC, dx$SCJ+dx$SCE), rownames(x))
        text(dx$IGC, dx$SCJ+dx$SCE, as.character(seq_len(nrow(dx))),
            cex=0.5, pos=pos)
    }
    ## IGC x SCE
    plot(dx$IGC, dx$SCE, panel.first=grid(), las=1, xlab=xlab,
        ylab="Splice-comp. exonic aln. (%)",col=grpcol[as.integer(group)], ...)
    if (labelpoints) {
        pos <- setNames(thigmophobe(dx$IGC, dx$SCJ+dx$SCE), rownames(x))
        text(dx$IGC, dx$SCJ+dx$SCE, as.character(seq_len(nrow(dx))),
            cex=0.5, pos=pos)
    }
    ## IGC x INT
    plot(dx$IGC, dx$INT, panel.first=grid(), las=1, xlab=xlab,
        ylab="Intronic alignments (%)", col=grpcol[as.integer(group)], ...)
    if (labelpoints) {
        pos <- setNames(thigmophobe(dx$IGC, dx$INT), rownames(x))
        text(dx$IGC, dx$INT, as.character(seq_len(nrow(dx))), cex=0.5, pos=pos)
    }
    ## STRANDEDNESS
    if (!all(is.na(dx$STRAND))) {
        plot(dx$IGC, dx$STRAND, panel.first=grid(), las=1, xlab=xlab,
            ylab="Strandedness", col=grpcol[as.integer(group)], ...)
        if (labelpoints) {
            pos <- setNames(thigmophobe(dx$IGC, dx$STRAND), rownames(x))
            text(dx$IGC, dx$STRAND, as.character(seq_len(nrow(dx))),
                cex=0.5, pos=pos)
        }
    }
    if (is.factor(group)) {
        par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
        plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
        legend("bottom", levels(group), fill=grpcol, inset=0,
            bg="white", xpd=NA, horiz = TRUE)
    }
})

## private function .setColorGrouping()

#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
.setColorGrouping <- function(group, dx) {
    grpcol <- "black"
    if (length(group) > 1) {
        if (is.numeric(group))
            group <- as.character(group)
        if (is.character(group) || is.factor(group)) {
            if (length(group) != nrow(dx)) {
                messgr<- paste("argument 'group' is a string character vector",
                        "or a factor, but it has a length different to the",
                        "number of BAM files in the input 'gDNAx' object.",
                        sep=" ")
                stop(messgr)
            }
            group <- factor(group)
            set1pal <- brewer.pal(nlevels(group), "Set1")
            if (nlevels(group) >= 6) {
                set1pal <- brewer.pal(nlevels(group)+1, "Set1")
                set1pal <- set1pal[-6] ## remove yellow, difficult to see
            }
            grpcol <- colorRampPalette(set1pal)(nlevels(group))

        } else
            stop("argument 'group' should be a string character vector or",
                "a factor, of the same length than the number of BAM files",
                "in the input 'gDNAx' object.")
    } else if (group != 1)
        stop("If 'group' has only one value, this value must be 1.")

    list(group=group, col=grpcol)
}

#' Plot alignment origins
#'
#' Using the output from gDNAdx(), plot the genomic origin of the alignments.
#'
#' @param x A 'gDNAx' object.
#'
#' @param group A string character vector or a factor, with as many values
#' as BAM files analyzed in 'x', whose values define groups among those
#' BAM files.
#'
#' @importFrom graphics barplot mtext par legend lines axis boxplot points
#' @importFrom stats density setNames
#' 
#' @examples
#' # plot origin of alignments per sample
#' plotAlnOrigins(gdnax, group=pdat$gDNA)
#' 
#' @export
#' @rdname gDNAdx
plotAlnOrigins <- function(x, group=1L) {
    if (nrow(getDx(x)) == 0)
        stop("no diagnostics available in the input 'gDNAx' object.")

    dx <- getDx(x)
    grp <- .setColorGrouping(group, dx)
    group <- grp$group
    grpcol <- grp$col
    if (length(group) > 1) {
        dx <- dx[order(group), ]
        group <- sort(group) ## group together samples from the same group
    }

    pct <- as.matrix(dx[, c("SCJ", "SCE", "INT", "IGC")])
    pct <- cbind(pct, OTH=100-rowSums(pct))
    alntypecolors <- c(IGC="darkblue", INT="skyblue",
                        SCJ="darkgreen", SCE="darkolivegreen",
                        OTH="darkgrey")
    bp <- barplot(t(pct), las=1, xaxt="n", col=alntypecolors[colnames(pct)],
                    ylab="Alignments (%)")
    mtext(rownames(pct), side=1, at=bp, las=2, cex=0.8,
            col=grpcol[as.integer(group)])
    par(xpd=TRUE)
    legend("topright", colnames(pct), bg="white",
            fill=alntypecolors[colnames(pct)], inset=c(-0.03, -0.03))
    par(xpd=FALSE)
}

#' Plot fragments length distributions
#'
#' Plot fragments length distributions estimated with gDNAdx()
#'
#' @param x A 'gDNAx' object.
#'
#' @examples
#' # plot fragments length distributions
#' plotFrgLength(gdnax)
#' 
#' @export
#' @rdname gDNAdx
plotFrgLength <- function(x) {
    alligcfrglen <- unlist(x@igcfrglen, use.names=FALSE)
    allintfrglen <- unlist(x@intfrglen, use.names=FALSE)
    allscjfrglen <- unlist(x@scjfrglen, use.names=FALSE)
    allscefrglen <- unlist(x@scefrglen, use.names=FALSE)
    igcden <- density(log10(alligcfrglen))
    intden <- density(log10(allintfrglen))
    scjden <- density(log10(allscjfrglen))
    sceden <- density(log10(allscefrglen))
    xrng <- range(c(igcden$x, intden$x, scjden$x, sceden$x))
    xrng[1] <- floor(xrng[1])
    xrng[2] <- floor(xrng[2])
    yrng <- range(c(igcden$y, intden$y, scjden$y, sceden$y))
    par(mfrow=c(1, 2), mar=c(4, 5, 1, 2))
    plot(igcden, lwd=2, col="darkblue", xaxt="n", xlab="log10 fragment length",
        las=1, main="", xlim=xrng, ylim=yrng)
    lines(intden, lwd=2, col="skyblue")
    lines(scjden, lwd=2, col="darkgreen")
    lines(sceden, lwd=2, col="darkolivegreen")
    ticks <- seq(xrng[1], xrng[2], by=1)
    axis(1, at=ticks, labels=parse(text=paste0("10^", ticks)))
    legend("topright", c("IGC", "INT", "SCJ", "SCE"), lwd=2, inset=0.01,
            col=c("darkblue", "skyblue", "darkgreen", "darkolivegreen"))
    yrng <- range(log10(c(alligcfrglen, allintfrglen, allscjfrglen,
                            allscefrglen)))
    yrng[1] <- floor(yrng[1])
    yrng[2] <- floor(yrng[2])
    allfrglen <- list(IGC=log10(alligcfrglen), INT=log10(allintfrglen),
                        SCJ=log10(allscjfrglen), SCE=log10(allscefrglen))
    boxplot(allfrglen, main="", yaxt="n", ylim=yrng,
            ylab="log10 fragment length",
            col=c("darkblue", "skyblue", "darkgreen", "darkolivegreen"))
    points((seq_len(length(allfrglen)))+0.25,
            c(mean(log10(alligcfrglen)), mean(log10(allintfrglen)),
            mean(log10(allscjfrglen)), mean(log10(allscefrglen))),
            pch=23, bg="black")
    ticks <- seq(yrng[1], yrng[2], by=1)
    axis(2, at=ticks, labels=parse(text=paste0("10^", ticks)), las=1)
}

