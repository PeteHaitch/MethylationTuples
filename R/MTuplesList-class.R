### =========================================================================
### GTuplesList objects
### -------------------------------------------------------------------------
###

#' MTuplesList objects
#' 
#' @description
#' The \code{MTuplesList} class is a container for the genomic locations 
#' of methylation loci in a genome.
#' 
#' @details
#' The \code{MTuplesList} class extends the \linkS4class{GTuplesList} class by 
#' adding the \code{methinfo} slot, which records the information about the 
#' type of methylation loci that are stored in the MTuplesList object.
#' 
#' @seealso 
#' \itemize{
#'  \item \linkS4class{GTuplesList} for the class from which MTuplesList 
#'  inherits
#'  \item \linkS4class{MTuples}
#' }
#' 
#' @aliases MTuplesList
#'
#' @include MethInfo-class.R
#' @author Peter Hickey
#' @examples
#' ## TODO
#' 
#' @importClassesFrom GenomicTuples GTuplesList
#' @importFrom methods setClass
#' 
#' @export
setClass("MTuplesList",
         contains = c("GTuplesList"),
         representation(
           unlistData = "MTuples",
           elementMetadata = "DataFrame"
         ),
         prototype(
           elementType = "MTuples"
         )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @rdname MTuplesList-class
#' @param ... \linkS4class{MTuples} objects
#' @importFrom IRanges PartitioningByEnd
#' @export
MTuplesList <- function(...) {
  listData <- list(...)
  if (length(listData) == 0L) {
    unlistData <- MTuples()
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]])) {
      listData <- listData[[1L]]
    }
    if (!all(sapply(listData, is, "MTuples"))) {
      stop("all elements in '...' must be MTuples objects")
    }
    if (!GenomicTuples:::.zero_range(sapply(listData, size)) && 
          !isTRUE(all(is.na(sapply(listData, size))))) {
      stop("all MTuples in '...' must have the same 'size'")
    }
    # TODO: Why suppressWarnings()?
    unlistData <- suppressWarnings(do.call("c", unname(listData)))
  }
  
  relist(unlistData, PartitioningByEnd(listData))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @rdname MTuplesList-class
#' @param x,object An MTuplesList object.
#' @importFrom methods setMethod
#' @export
setMethod("methinfo", "MTuplesList", 
          function(x) {
            x@unlistData@methinfo
          }
)

#' @rdname MTuplesList-class
#' @inheritParams methinfo
#' @importFrom methods setMethod
#' @export
setMethod("methtype", "MTuplesList", 
          function(x) {
            methtype(x@unlistData@methinfo)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

#' @rdname MTuplesList-class
#' @inheritParams methinfo
#' @param value A replacement object of the appropriate class.
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methinfo", c("MTuplesList", "MethInfo"), 
                 function(x, value) {
                   x@unlistData@methinfo <- value
                   x
                 }
)

#' @rdname MTuplesList-class
#' @inheritParams methinfo<-
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methtype", c("MTuplesList", "character"), 
                 function(x, value) {
                   methtype(x@unlistData@methinfo) <- value
                   x
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from MTuples to MTuplesList with extractList() and family.
###

#' @rdname MTuplesList-class
#' @inheritParams methinfo
#' @importFrom methods setMethod
#' @importMethodsFrom IRanges relistToClass
#' @export
setMethod("relistToClass", "MTuples", 
          function(x) {
            "MTuplesList"
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

# Based on GenomicRanges::showList
#' @importFrom utils head tail
#' @importMethodsFrom GenomeInfoDb seqinfo
#' @importMethodsFrom S4Vectors elementLengths
my_showList <- function(object, showFunction, print.classinfo) {
  k <- length(object)
  cumsumN <- cumsum(elementLengths(object))
  N <- tail(cumsumN, 1)
  cat(class(object), " object of length ", k, ":\n", sep = "")
  if (k == 0L) {
    cat("<0 elements>\n\n")
  } else if ((k == 1L) || ((k <= 3L) && (N <= 20L))) {
    nms <- names(object)
    defnms <- paste0("[[", seq_len(k), "]]")
    if (is.null(nms)) {
      nms <- defnms
    } else {
      empty <- nchar(nms) == 0L
      nms[empty] <- defnms[empty]
      nms[!empty] <- paste0("$", nms[!empty])
    }
    for (i in seq_len(k)) {
      cat(nms[i], "\n")
      showFunction(object[[i]], margin = "  ", 
                   print.classinfo = print.classinfo)
      if (print.classinfo) {
        print.classinfo <- FALSE
      }
      cat("\n")
    }
  } else {
    sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
    if (k >= 3 && cumsumN[3L] <= 20) {
      showK <- 3
    }
    else if (k >= 2 && cumsumN[2L] <= 20) {
      showK <- 2
    }
    else {
      showK <- 1
    }
    diffK <- k - showK
    nms <- names(object)[seq_len(showK)]
    defnms <- paste0("[[", seq_len(showK), "]]")
    if (is.null(nms)) {
      nms <- defnms
    } else {
      empty <- nchar(nms) == 0L
      nms[empty] <- defnms[empty]
      nms[!empty] <- paste0("$", nms[!empty])
    }
    for (i in seq_len(showK)) {
      cat(nms[i], "\n")
      showFunction(object[[i]], margin = "  ",
                   print.classinfo = print.classinfo)
      if (print.classinfo) {
        print.classinfo <- FALSE
      }
      cat("\n")
    }
    if (diffK > 0) {
      cat("...\n<", k - showK,
          ifelse(diffK == 1, " more element>\n", " more elements>\n"), sep = "")
    }
  }
  cat("-------\n")
  cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
  cat("methinfo: ", summary(methinfo(object)), "\n", sep = "")
}

#' @rdname MTuplesList-class
#' @inheritParams methinfo
#' @importFrom methods setMethod
#' @export
setMethod("show", "MTuplesList",
          function(object) {
            my_showList(object, showMTuples, FALSE)
          }
)
