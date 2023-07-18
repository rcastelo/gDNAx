library(gDNAinRNAseqData)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


test_input_errors_identifyStrandMode <- function() {

    # Retrieving BAM files
    bamfiles <- LiYu22subsetBAMfiles()
    
    ## bamfiles
    bamfiles2 <- c("invented_sample.bam")
    checkException(identifyStrandMode(bamfiles2, txdb, singleEnd=FALSE),
                    "testing bamfiles input error")
    
    ## txdb
    txdb2 <- "invented_txdb_packagename"
    checkException(identifyStrandMode(bamfiles, txdb2, singleEnd=FALSE),
                    "testing txdb input error")
    
    ## singleEnd
    checkException(identifyStrandMode(bamfiles, txdb, singleEnd="no"),
                    "testing singleEnd input error")
    
    ## stdChrom
    checkException(identifyStrandMode(bamfiles, txdb, singleEnd=FALSE,
                        stdChrom="yes"), "testing stdChrom input error")
    
    ## yieldSize
    checkException(identifyStrandMode(bamfiles, txdb, singleEnd=FALSE,
                        yieldSize = "100"), "testing yieldSize input error")
    checkException(identifyStrandMode(bamfiles, txdb, singleEnd=FALSE,
                        yieldSize = NA), "testing yieldSize input error")
    checkException(identifyStrandMode(bamfiles, txdb, singleEnd=FALSE,
                        yieldSize = 100), "testing yieldSize input error")
}

test_output_identifyStrandMode <- function() {

    # Retrieving BAM files
    bamfiles <- LiYu22subsetBAMfiles()
    
    suppressWarnings(strandM <- identifyStrandMode(bamfiles, txdb,
                                                    singleEnd=FALSE))
    
    checkTrue(is(strandM, "list"))
    checkEquals(length(strandM), 2)
    checkTrue(is(strandM$Strandedness, "matrix"))
    
    # All % and frag length values should be numeric and between 0 and 1
    checkTrue(is(unlist(as.vector(strandM$Strandedness[,1:3])), "numeric"))
    checkTrue(all(strandM$Strandedness[,1:3] >= 0 & 
                    strandM$Strandedness[,1:3] <= 1))
    
    # strandMode1 + strandMode2 + ambig should add 1
    checkEquals(as.vector(rowSums(strandM$Strandedness[,1:3])),
                rep(1, length(bamfiles)), tolerance = 1.0e-3)
}

test_decideStrandMode <- function() {

    ## ambiguous
    strbysm <- data.frame("strandMode1" = 0.3, "strandMode2" = 0.65,
                            "ambig" = 0.05)
    checkIdentical(gDNAx:::.decideStrandMode(strbysm), "ambiguous")
    
    ## strandMode2
    strbysm <- data.frame("strandMode1" = 0.01, "strandMode2" = 0.94,
                          "ambig" = 0.05)
    checkIdentical(gDNAx:::.decideStrandMode(strbysm), 2L)
    
    ## strandMode1
    strbysm <- data.frame("strandMode1" = 0.91, "strandMode2" = 0.08,
                          "ambig" = 0.01)
    checkIdentical(gDNAx:::.decideStrandMode(strbysm), 1L)
    
    ## unstranded: NA
    strbysm <- data.frame("strandMode1" = 0.55, "strandMode2" = 0.42,
                          "ambig" = 0.02)
    checkTrue(is.na(gDNAx:::.decideStrandMode(strbysm)))
}
