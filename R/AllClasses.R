#' gDNAx class
#'
#' This is a class for storing the results of a call to the
#' 'gDNAdx()' function.
#'
#' @slot bfl A \linkS4class{BamFileList} object.
#'
#' @slot txdbpkg A \linkS4class{TxDb} object.
#'
#' @slot singleEnd Logical value indicating if reads are single (\code{TRUE})
#' or paired-end (\code{FALSE}).
#'
#' @slot strandMode Numeric vector which can take values 0, 1 or 2. The strand
#' mode is a per-object switch on
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' objects that controls the behavior of the strand getter. See
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail.
#'
#' @slot stdChrom Logical value indicating whether only alignments in the
#' 'standard chromosomes' should be used. Consult the help page of the function
#' \code{\link[GenomeInfoDb]{keepStandardChromosomes}} from the package
#' \code{GenomeInfoDb} for further information.
#'
#' @slot readLength Integer value storing the read length.
#'
#' @slot yieldSize Integer value storing the number of alignments employed by
#' the function \code{\link{gDNAdx}()}.
#'
#' @slot diagnostics A 'data.frame' object storing the diagnostics calculated
#' by the function 'gDNAdx()'.
#'
#' @slot igcfrglen A 'list' object storing the fragment lengths derived from
#' alignments in intergenic regions.
#'
#' @slot intfrglen A 'list' object storing the fragment lengths derived from
#' alignments in intronic regions.
#'
#' @slot scjfrglen A 'list' object storing the fragment lengths derived from
#' spliced-compatible junction alignments in transcripts.
#'
#' @slot scefrglen A 'list' object storing the fragment lengths derived
#' from spliced-compatible exonic alignments in transcripts.
#'
#' @slot sicfrglen A 'list' object storing the fragment lengths derived
#' from splice-incompatible alignments in transcripts.
#'
#' @slot intergenic A 'GRanges' object storing the intergenic feature
#' annotations.
#'
#' @slot intronic A 'GRanges' object storing the intronic feature
#' annotations.
#'
#' @slot transcripts A 'GRangesList' object storing the transcript annotations.
#'
#' @slot tx2gene A string character vector storing the correspondence between
#' transcripts and genes according to an 'TxDb' object.
#'
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Retrieving BAM files
#' bamfiles <- LiYu22subsetBAMfiles()
#' 
#' gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' gdnax
#'
#' @name gDNAx-class
#' @rdname gDNAx-class
#' @exportClass gDNAx
setClass("gDNAx",
        representation(bfl="BamFileList",
                        txdbpkg="character",
                        singleEnd="logical",
                        strandMode="integer",
                        stdChrom="logical",
                        readLength="integer",
                        yieldSize="integer",
                        diagnostics="data.frame",
                        igcfrglen="list",
                        intfrglen="list",
                        scjfrglen="list",
                        scefrglen="list",
                        sicfrglen="list",
                        intergenic="GRanges",
                        intronic="GRanges",
                        transcripts="GRangesList",
                        tx2gene="character"))

#' @param x A \linkS4class{gDNAx} object.
#'
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Getting the 'gDNAx' object. Can be done using the commented code
#' # bamfiles <- LiYu22subsetBAMfiles()
#' # gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' 
#' # Here to reduce example running time, the 'gDNAx' object is loaded
#' gdnax <- file.path(system.file("extdata", package="gDNAx"), "gdnax.rds")
#' 
#' # Getting statistics
#' dx <- getDx(gdnax)
#' head(dx)
#'
#' @export
#' @aliases getDx
#' @aliases getDx,gDNAx-method
#' @rdname gDNAx-class
#' @name getDx
setMethod("getDx", "gDNAx",
            function(x) {
                x@diagnostics
            })

#' @param object A \linkS4class{gDNAx} object.
#'
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Getting the 'gDNAx' object. Can be done using the commented code
#' # bamfiles <- LiYu22subsetBAMfiles()
#' # gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' 
#' # Here to reduce example running time, the 'gDNAx' object is loaded
#' gdnax <- file.path(system.file("extdata", package="gDNAx"), "gdnax.rds")
#' gdnax
#' 
#' @importFrom methods show
#' 
#' @export
#' @aliases show
#' @aliases show,gDNAx-method
#' @rdname gDNAx-class
setMethod("show", "gDNAx",
            function(object) {
                cat(class(object), "object\n")
                cat(sprintf("# BAM files (%d): %s\n", length(object@bfl),
                            .pprintnames(names(object@bfl))))
                if (object@singleEnd)
                    cat(sprintf("# Library layout: single-end (%dnt)\n",
                                object@readLength))
                else {
                    cat(sprintf("# Library layout: paired-end (2x%dnt)\n",
                                object@readLength))
                    cat(sprintf("# Strand mode: %d\n", object@strandMode))
                }
                if (object@stdChrom)
                    cat("# Sequences: only standard chromosomes\n")
                else
                    cat("# Sequences: all\n")
                cat(sprintf("# Annotation pkg: %s\n", object@txdbpkg))
                cat(sprintf("# Alignments employed: first %d\n",
                            object@yieldSize))
            })

#' @param x A \linkS4class{gDNAx} object.
#'
#' @return \code{features()}: A \code{GRanges} object with intergenic ranges.
#' 
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Getting the 'gDNAx' object. Can be done using the commented code
#' # bamfiles <- LiYu22subsetBAMfiles()
#' # gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' 
#' # Here to reduce example running time, the 'gDNAx' object is loaded
#' gdnax <- file.path(system.file("extdata", package="gDNAx"), "gdnax.rds")
#' 
#' igc <- getIgc(gdnax)
#' head(igc, n=3)
#' 
#' @export
#' @aliases getIgc
#' @aliases getIgc,gDNAx-method
#' @rdname gDNAx-class
#' @name getIgc
setMethod("getIgc", "gDNAx",
            function(x) {
                x@intergenic
            })

#' @param x A \linkS4class{gDNAx} object.
#'
#' @return \code{features()}: A \code{GRanges} object with intron ranges.
#' 
#' @examples
#' library(gDNAinRNAseqData)
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' 
#' # Getting the 'gDNAx' object. Can be done using the commented code
#' # bamfiles <- LiYu22subsetBAMfiles()
#' # gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA)
#' 
#' # Here to reduce example running time, the 'gDNAx' object is loaded
#' gdnax <- file.path(system.file("extdata", package="gDNAx"), "gdnax.rds")
#' 
#' int <- getInt(gdnax)
#' head(int, n=3)
#' 
#' @export
#' @aliases getInt
#' @aliases getInt,gDNAx-method
#' @rdname gDNAx-class
#' @name getInt
setMethod("getInt", "gDNAx",
            function(x) {
                x@intronic
            })
