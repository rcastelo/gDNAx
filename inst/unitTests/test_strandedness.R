library(gDNAinRNAseqData)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Retrieving BAM files
bamfiles <- LiYu22subsetBAMfiles()

# Getting information about the gDNA concentrations of each BAM file
pdat <- LiYu22phenoData(bamfiles)


test_input_errors_strandedness <- function() {

    # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    
    ## bamfiles
    bamfiles2 <- c("invented_sample.bam")
    checkException(strandedness(bamfiles2, txdb, singleEnd=FALSE),
                    "testing bamfiles input error")
    
    ## txdb
    txdb2 <- "invented_txdb_packagename"
    checkException(strandedness(bamfiles, txdb2, singleEnd=FALSE),
                    "testing txdb input error")
    
    ## singleEnd
    checkException(strandedness(bamfiles, txdb, singleEnd="no"),
                    "testing singleEnd input error")
    
    ## stdChrom
    checkException(strandedness(bamfiles, txdb, singleEnd=FALSE, stdChrom="yes"),
                   "testing stdChrom input error")
    
}

test_output_strandedness <- function() {

    # Retrieving BAM files
    # bamfiles <- LiYu22subsetBAMfiles()
    
    suppressWarnings(strness <- strandedness(bamfiles[1:2], txdb, singleEnd=FALSE))
    
    checkTrue(is(strness, "data.frame"))
    
    # All % and frag length values should be numeric and between 0 and 1
    checkTrue(is(unlist(strness[, 1:3]), "numeric"))
    checkTrue(all(unlist(strness[, 1:3]) >= 0 & unlist(strness[, 1:3]) <= 1))
    
    # strandMode1 + strandMode2 + ambig should add 1
    checkEquals(as.vector(rowSums(strness[, 1:3])),
                rep(1, length(bamfiles[1:2])), tolerance = 1.0e-3)
}

test_classifyStrandMode <- function() {

    ## unstranded: NA
    strbysm <- data.frame("strandMode1" = 0.3, "strandMode2" = 0.45,
                          "ambig" = 0.25)
    checkTrue(is.na(classifyStrandMode(strbysm)))
    
    ## strandMode2
    strbysm <- data.frame("strandMode1" = 0.01, "strandMode2" = 0.94,
                          "ambig" = 0.05)
    checkIdentical(as.integer(classifyStrandMode(strbysm)), 2L)
    
    ## strandMode1
    strbysm <- data.frame("strandMode1" = 0.91, "strandMode2" = 0.08,
                          "ambig" = 0.01)
    checkIdentical(as.integer(classifyStrandMode(strbysm)), 1L)
    
}
