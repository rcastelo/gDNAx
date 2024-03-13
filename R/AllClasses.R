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
#' @slot strandMode Integer value either 0, 1, 2 or \code{NA}, indicating how
#' the strand of a pair of read alignments should be inferred from the strand
#' of the first and last alignments in the pair. See the
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail.
#'
#' @slot allStrandModes Vector of integer values each of them corresponding to
#' a \code{strandMode} value estimated from a BAM file.
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
#' @slot strandedness A 'data.frame' object storing the estimated values of
#' strandedness, calculated when the argument \code{strandMode} is missing in
#' the call to the function 'gDNAdx()'.
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
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' bamfiles <- LiYu22subsetBAMfiles() # Retrieving BAM files
#' ## one could simply call 'gDNAx(bamfiles, txdb)' but we give the arguments
#' ## below to reduce time and verbosity when running this example
#' gdnax <- gDNAdx(bamfiles, txdb, singleEnd=FALSE, strandMode=NA,
#'                 useRMSK=FALSE, verbose=FALSE)
#' gdnax
#' 
#'
#' @name gDNAx-class
#' @rdname gDNAx-class
#' @exportClass gDNAx
setClass("gDNAx",
         representation(bfl="BamFileList",
                        txdbpkg="character",
                        singleEnd="logical",
                        strandMode="integer",
                        allStrandModes="integer",
                        stdChrom="logical",
                        readLength="integer",
                        yieldSize="integer",
                        diagnostics="data.frame",
                        strandedness="data.frame",
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
                rlen <- object@readLength
                urlen <- sort(unique(rlen))
                rlenstr <- sprintf("%d (%s%dnt)", table(rlen),
                                   ifelse(object@singleEnd, "", "2x"),
                                   urlen)
                rlenstr <- sprintf("%s %s",
                                   ifelse(object@singleEnd, "single-end,",
                                          "paired-end,"),
                                   paste(rlenstr, collapse=", "))
                cat(sprintf("# Library layout: %s\n", rlenstr))
                pstr <- smstr <- ""
                if (length(object@allStrandModes) > 1) {
                    tab <- table(object@allStrandModes, useNA="always")
                    maxsm <- names(which.max(tab))
                    nmaxsm <- tab[maxsm]
                    if (is.na(maxsm)) ## unstranded
                      nmaxsm <- sum(is.na(object@allStrandModes))
                    pstr <- sprintf("(%d out of %d)", nmaxsm, sum(tab))
                    if (!is.na(maxsm)) ## stranded
                        pstr <- sprintf("(%d out of %d)",
                                        sum(tab[!is.na(names(tab))]), sum(tab))
                    smstr <- sprintf("(%d in mode 1, %d in mode 2)", tab["1"],
                                     tab["2"])
                }
                if (is.na(object@strandMode))
                    cat(sprintf("# Library protocol: unstranded %s\n", pstr))
                else {
                    if (object@singleEnd)
                        cat(sprintf("# Library protocol: stranded %s\n", pstr))
                    else {
                        ststr <- paste("# Library protocol: stranded %s,",
                                       "mode %d %s\n")
                        cat(sprintf(ststr, pstr, object@strandMode, smstr))
                    }
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

#' @param x A \linkS4class{gDNAx} object.
#'
#' @examples
#' singleEnd(gdnax)
#' 
#' @export
#' @aliases singleEnd
#' @aliases singleEnd,gDNAx-method
#' @rdname gDNAx-class
setMethod("singleEnd", "gDNAx",
          function(x) x@singleEnd
         )

#' @param x A \linkS4class{gDNAx} object.
#'
#' @examples
#' strandMode(gdnax)
#' 
#' @importFrom GenomicAlignments strandMode
#' 
#' @export
#' @aliases strandMode
#' @aliases strandMode,gDNAx-method
#' @rdname gDNAx-class
setMethod("strandMode", "gDNAx",
          function(x) x@strandMode
         )

#' @param x A \linkS4class{gDNAx} object.
#' @param value Integer value either 0, 1, 2 or \code{NA}, indicating how
#' the strand of a pair of read alignments should be inferred from the strand
#' of the first and last alignments in the pair. See the
#' \code{\link[GenomicAlignments:GAlignmentPairs-class]{GAlignmentPairs}}
#' class for further detail.
#'
#' @examples
#' strandMode(gdnax) <- NA
#' 
#' @importFrom S4Vectors isSingleNumber
#' @importFrom GenomicAlignments strandMode<-
#' 
#' @export
#' @aliases strandMode<-
#' @aliases strandMode<-,gDNAx-method
#' @rdname gDNAx-class
setReplaceMethod("strandMode", "gDNAx",
                 function(x, value) {
                     if (singleEnd(x))
                         stop("Cannot set strand mode on single-end data.")
                     if (!is.na(value)) {
                         if (!isSingleNumber(value))
                           stop("invalid strand mode (must be 0, 1, 2, or NA).")
                         if (!is.integer(value))
                           value <- as.integer(value)
                         if (!(value %in% 0:2))
                           stop("invalid strand mode (must be 0, 1, 2, or NA).")
                     } else
                         value <- as.integer(value)
                     x@strandMode <- value
                     x
                 })

#' @param x A \linkS4class{gDNAx} object.
#'
#' @examples
#' allStrandModes(gdnax)
#' 
#' @export
#' @aliases allStrandModes
#' @aliases allStrandModes,gDNAx-method
#' @rdname gDNAx-class
setMethod("allStrandModes", "gDNAx",
          function(x) x@allStrandModes
         )

#' @param x A \linkS4class{gDNAx} object.
#'
#' @examples
#' strandedness(gdnax)
#' 
#' @export
#' @aliases strandedness
#' @aliases strandedness,gDNAx-method
#' @rdname gDNAx-class
setMethod("strandedness", "gDNAx",
          function(x) x@strandedness
         )

