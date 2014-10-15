### =========================================================================
### MTuples objects: genomic tuples of methylation loci
### -------------------------------------------------------------------------
###

# TODO: Methods to create MTuples from GRanges and BSgenome objects.
#' MTuples objects
#' 
#' @description
#' The \code{MTuples} class is a container for the genomic locations 
#' of methylation loci in a genome.
#' 
#' @details
#' The \code{MTuples} class extends the 
#' \code{\link[GenomicTuples]{GTuples}} class by adding the 
#' \code{methinfo} slot, which records the information about the type of 
#' methylation loci that are stored in the object. 
#' 
#' @param methinfo A \code{\link{MethInfo}} object containing information about 
#' the methylation loci present in the \code{MTuples} object.
#' @param gtuples A \code{\link{GTuples}} object containing the positions of 
#' the methylation loci as genomic tuples.
#' 
#' @seealso \code{\link{findMTuples}} to find genomic tuples of methylation 
#' loci in a \code{\link[BSgenome]{BSgenome}} object, and
#' \code{\link[GenomicTuples]{GTuples}} for the class from which \code{MTuples} 
#' inherits.
#'
#' @export
#' @include MethInfo-class.R
setClass("MTuples",
         contains = "GTuples",
         representation(
           methinfo = "MethInfo"),
         prototype(
           methinfo = MethInfo())
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

# None as yet. Validity methods "inherited" from GRanges and MethInfo classes.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
MTuples <- function(gtuples = GTuples(), methinfo = MethInfo()) {
  new("MTuples", gtuples, methinfo = methinfo)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

# TODO: Need to keep this up to date with the show,GTuples-method
showMTuples <- function(x, margin = "", print.classinfo = FALSE, 
                        print.seqinfo = TRUE, print.methinfo = TRUE) {
  if (!identical(print.classinfo, FALSE)) {
    stop("'print.classinfo' not implemented")
  }
  lx <- length(x)
  nc <- ncol(mcols(x))
  
  if (!is.na(x@size)) {
    cat(class(x), " object with ", lx, " x ", 
        ifelse(lx == 1L, paste0(x@size, "-tuple"), 
               paste0(x@size, "-tuples")), 
        " and ", nc, " metadata ", ifelse(nc == 1L, "column", "columns"), 
        ":\n", sep = "")
  } else{
    cat(class(x), " with 0 tuples and 0 metadata columns:\n", sep = "")
  }
  
  out <- S4Vectors:::makePrettyMatrixForCompactPrinting(
    x, GenomicTuples:::.makeNakedMatFromGTuples)
  # TODO: Try to implement 'print.classinfo', although low priority.
  ## These lines commented out because classinfo is more complicated for GTuples 
  ## objects than GRanges objects. For example, some of the `pos` information 
  ## is stored in an IRanges object while some is stored in a matrix.
  #if (print.classinfo) {
  #    .COL2CLASS <- c(seqnames = "Rle", ranges = "IRanges", 
  #        strand = "Rle")
  #    extraColumnNames <- extraColumnSlotNames(x)
  #    .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
  #    classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
  #    stopifnot(identical(colnames(classinfo), colnames(out)))
  #    out <- rbind(classinfo, out)
  #}
  
  if (nrow(out) != 0L){ 
    rownames(out) <- paste0(margin, rownames(out))
  }
  print(out, quote = FALSE, right = TRUE)
  if (print.seqinfo & print.methinfo) {
    cat(margin, "---\n", sep = "")
    cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep ="")
    cat(margin, "methinfo: ", summary(methinfo(x)), "\n", sep = "")
  } else if (print.seqinfo) {
    cat(margin, "---\n", sep = "")
    cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep ="")
  } else if (print.methinfo) {
    cat(margin, "---\n", sep = "")
    cat(margin, "methinfo: ", summary(methinfo(x)), "\n", sep = "")
  }
}

#' @export
setMethod("show", 
          "MTuples", 
          function(object) {
            showMTuples(object, margin="  ", print.seqinfo = TRUE, 
                        print.methinfo = TRUE)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

.unlist_list_of_MTuples <- function(x, ignore.mcols = FALSE) {
  if (!isTRUEorFALSE(ignore.mcols)) {
    stop("'ignore.mcols' must be TRUE or FALSE")
  }
  ans_class <- class(x[[1L]])
  ans_methinfo <- do.call(merge, lapply(x, methinfo))
  ans_seqinfo <- do.call(merge, lapply(x, seqinfo))
  ans_seqnames <- do.call(c, lapply(x, seqnames))
  ans_tuples <- unname(do.call(rbind, lapply(x, tuples)))
  ans_strand <- do.call(c, lapply(x, strand))
  if (ignore.mcols) {
    ans_mcols <- new("DataFrame", nrows = length(ans_tuples))
  } else {
    ans_mcols <- do.call(rbind, lapply(x, mcols, FALSE))
  }
  new(ans_class, 
      GTuples(seqnames = ans_seqnames, tuples = ans_tuples, 
      strand = ans_strand, ans_mcols, seqinfo = ans_seqinfo),
      methinfo = ans_methinfo)
}

#' @export
setMethod("c", 
          "MTuples", 
          function(x, ..., ignore.mcols = FALSE, recursive = FALSE) {
            if (!identical(recursive, FALSE)) {
              stop("'recursive' argument not supported")
            }
            if (missing(x)) {
              args <- unname(list(...))
            } else {
              args <- unname(list(x, ...))
            }
            if (!GenomicTuples:::.zero_range(sapply(args, size)) && 
                  !isTRUE(all(is.na(sapply(args, size))))) {
              stop(paste0("Cannot combine ", paste0(unique(sapply(args, class)), 
                                                    collapse = ' and '), 
                          " containing tuples of different 'size'."))
            }
            .unlist_list_of_MTuples(args, ignore.mcols = ignore.mcols)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

# TODO: methinfo, methtype
#' @export
setMethod("methinfo", 
          "MTuples", 
          function(x) {
            x@methinfo
          }
)

#' @export
setMethod("methtype", 
          "MTuples", 
          function(x) {
            methtype(x@methinfo)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

#' @export
setReplaceMethod("methinfo", 
                 c("MTuples", "MethInfo"), 
                 function(x, value) {
                   x@methinfo <- value
                   x
                 }
)

#' @export
setReplaceMethod("methtype", 
                 c("MTuples", "character"), 
                 function(x, value) {
                   methtype(x@methinfo) <- value
                   x
                 }
)
