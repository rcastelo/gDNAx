library(gDNAinRNAseqData)
library(BiocGenerics)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Retrieving BAM files
bamfiles <- LiYu22subsetBAMfiles()

# Getting information about the gDNA concentrations of each BAM file
pdat <- LiYu22phenoData(bamfiles)

test_download_LiYu22 <- function() {
    
    # # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    # 
    # # Getting information about the gDNA concentrations of each BAM file
    # pdat <- LiYu22phenoData(bamfiles)
    
    checkEquals(length(bamfiles), 9L)
    checkEquals(dim(pdat), c(9,1))
    checkEquals(length(unique(pdat[,1])), 3)
    
    bn <- gsub(".bam$", "", basename(bamfiles))
    checkIdentical(rownames(pdat), bn)
}

# test_output_gDNAdx <- function() {
#     
#     # # Retrieving BAM files
#     # bamfiles <- LiYu22subsetBAMfiles()
#     
#     gdnax <- gDNAdx(bamfiles[1], txdb, singleEnd=FALSE, strandMode=NA)
#     dx <- getDx(gdnax)
#     
#     checkTrue(is(gdnax, "gDNAx"))
#     checkEquals(dim(dx), c(1,10))
#     checkTrue(is(dx, "data.frame"))
#     
#     # If strandMode=NA, output of strandedness should be NA
#     checkTrue(is.na(unique(dx[,"STRAND"])))
#     
#     # All % and frag length values should be numeric and positive
#     checkTrue(is(unlist(as.vector(dx[,-10])), "numeric"))
#     checkTrue(all(dx[,-10] > 0))
#     
#     # igc + int + scj + sce should never add > 100%
#     checkTrue(all(rowSums(dx[,1:4]) < 100))
# }

test_input_errors <- function() {
    
    # # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    
    ## bamfiles
    bamfiles2 <- c("invented_sample.bam")
    checkException(gDNAdx(bamfiles2, txdb, singleEnd=FALSE, strandMode=NA),
                    "testing bamfiles input error")
    
    ## txdb
    txdb2 <- "invented_txdb_packagename"
    checkException(gDNAdx(bamfiles, txdb2, singleEnd=FALSE, strandMode=NA),
                    "testing txdb input error")
    
    ## singleEnd
    checkException(gDNAdx(bamfiles, txdb, singleEnd="no", strandMode=NA),
                    "testing singleEnd input error")
    
    ## stdChrom
    checkException(gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA,
                            stdChrom="yes"), "testing stdChrom input error")
    
    ## yieldSize
    checkException(gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA,
                        yieldSize = "100"), "testing yieldSize input error")
    checkException(gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA,
                        yieldSize = NA), "testing yieldSize input error")
    checkException(gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA,
                        yieldSize = 100), "testing yieldSize input error")
}

library(BiocParallel)
library(GenomicAlignments)
library(AnnotationHub)
library(GenomeInfoDb)
library(GenomicRanges)

test_minFrgLen <- function() {
    
    # # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    
    bfl <- gDNAx:::.checkBamFileListArgs(bamfiles[1:2], singleEnd=FALSE,
                                         fragments=FALSE, yieldSize=100000L)
    strandMode <- gDNAx:::.checkStrandMode(NA)
    minfrglen <- gDNAx:::.minFrgLen(bfl, singleEnd=FALSE)
    
    checkEquals(minfrglen, 100)
}

test_dx_computation <- function() {
    
    # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    
    bfl <- bamfiles[1]
    singleEnd=FALSE
    strandMode=NA
    stdChrom=TRUE
    yieldSize=100000L
    verbose=FALSE
    
    yieldSize <- gDNAx:::.checkYieldSize(yieldSize)
    bfl <- gDNAx:::.checkBamFileListArgs(bfl, singleEnd, fragments=FALSE,
                                            yieldSize)
    strandMode <- gDNAx:::.checkStrandMode(strandMode)
    gDNAx:::.checkOtherArgs(singleEnd, stdChrom)
    if (is.character(txdb))
        txdb <- gDNAx:::.loadAnnotationPackageObject(txdb, "txdb", "TxDb",
                                                verbose=verbose)
    if (verbose)
        message(sprintf("Fetching annotations for %s", genome(txdb)[1]))
    minfrglen <- gDNAx:::.minFrgLen(bfl, singleEnd)
    igcintrng <- gDNAx:::.fetchIGCandINTrng(txdb, minfrglen, stdChrom,
                                            strandMode)
    
    checkEquals(length(igcintrng), 2)
    checkIdentical(names(igcintrng), c("igcrng", "intrng"))
    
    exbytx <- exonsBy(txdb, by="tx") ## fetch transcript annotations
    if (stdChrom) {
        exbytx2 <- keepStandardChromosomes(exbytx, pruning.mode="fine")
        exbytx2 <- exbytx2[lengths(exbytx2) > 0]
    }
    
    ## Checking reduction on number of ranges after removing ranges in
    ## non-standard chromosomes
    checkTrue(length(exbytx) >= length(exbytx2))
    exbytx <- exbytx2
    
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
    
    ## Using only 1 file for testing
    bf <- bfl[[1]]
    igc=igcintrng$igcrng
    int=igcintrng$intrng
    tx=exbytx
    
    ## Running .gDNAdx_oneBAM() for this bam file
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
      gal <- keepStandardChromosomes(gal, pruning.mode="fine")
    
    if (verbose)
      message(sprintf("Processing alignments from %s", basename(path(bf))))
    
    gal <- gDNAx:::.matchSeqinfo(gal, tx, verbose)
    
    ## Testing that .matchSeqinfo() unified seqlevels
    checkIdentical(seqlevels(gal), seqlevels(tx))
    
    naln <- length(gal)
    
    ## intergenic alignments
    igcaln <- gDNAx:::.igcAlignments(gal, igc, fragmentsLen=TRUE, nfrgs=1000)
    ## intronic alignments
    intaln <- gDNAx:::.intAlignments(gal, int, strandMode, fragmentsLen=TRUE,
                                        nfrgs=1000)
    ## splice-compatible alignments
    scoaln <- gDNAx:::.scoAlignments(gal, tx, tx2gene, singleEnd, strandMode,
                                        fragmentsLen=TRUE, nfrgs=1000)
    
    ## Assigning reads simultaneously identified as IGC/INT and SCE/SCJ
    ## as SCE/SCJ
    whigc <- (igcaln$igcmask & scoaln$scjmask) | 
      (igcaln$igcmask & scoaln$scemask)
    igcaln$igcmask[whigc] <- FALSE
    
    whint <- (intaln$intmask & scoaln$scjmask) | 
      (intaln$intmask & scoaln$scemask)
    intaln$intmask[whint] <- FALSE
    
    ## Checking there are no reads simultaneously identified as IGC/INT and
    ## SCE/SCJ
    checkEquals(sum(intaln$intmask & scoaln$scjmask), 0)
    checkEquals(sum(intaln$intmask & scoaln$scemask), 0)
    checkEquals(sum(igcaln$igcmask & scoaln$scjmask), 0)
    checkEquals(sum(intaln$intmask & scoaln$scemask), 0)
    
    dx_oneBAM <- gDNAx:::.getdx_oneBAM(igcaln, intaln, scoaln, singleEnd,
                               strandMode, gal, tx, naln)
    
    checkEquals(length(dx_oneBAM), 11)
    checkEquals(dx_oneBAM$naln, 99660)
    checkEquals(dx_oneBAM$nigcaln, 1039)
    checkEquals(dx_oneBAM$nintaln, 11715)
    checkEquals(dx_oneBAM$nscjaln, 15129)
    checkEquals(dx_oneBAM$nscealn, 39942)
    checkEquals(dx_oneBAM$nsccaln, 19498)
}





