library(gDNAinRNAseqData)
library(gDNAx)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

bamfiles <- LiYu22subsetBAMfiles()
pdat <- LiYu22phenoData(bamfiles)

test_download_LiYu22 <- function() {
    
    checkEquals(length(bamfiles), 9L)
    checkEquals(dim(pdat), c(9L, 1L))
    checkEquals(length(unique(pdat[, 1])), 3L)
    
    bn <- gsub(".bam$", "", basename(bamfiles))
    checkIdentical(rownames(pdat), bn)
}

test_input_errors <- function() {
    
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
}

test_readLengths <- function() {
    
    bfl <- gDNAx:::.checkBamFileListArgs(bamfiles[1:2], singleEnd=FALSE, yieldSize=100000L)
    strandMode <- gDNAx:::.checkStrandMode(NA)
    rlens <- gDNAx:::.readLengths(bfl, singleEnd=FALSE, verbose=FALSE)
    
    checkEquals(unique(rlens), 50)
}

test_gDNAdx <- function() {

    gdnax <- gDNAdx(bamfiles[1:2], txdb, verbose=FALSE)

    checkEquals(unname(rowSums(strandedness(gdnax)[, 1:3])), c(1L, 1L))

    checkEquals(strandMode(gdnax), as.integer(NA))

    checkEquals(nrow(getDx(gdnax)), 2L)

     # igc + int + scj + sce should never add > 100%
    checkTrue(all(rowSums(getDx(gdnax)[, c("IGC", "INT", "SCJ", "SCE")]) < 100))
}
