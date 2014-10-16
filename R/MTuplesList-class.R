### =========================================================================
### GTuplesList objects
### -------------------------------------------------------------------------
###

# TODO: unit tests
# TODO: Base documentation on GTuplesList
#' MTuplesList objects
#' 
#' @description
#' The \code{MTuplesList} class is a container for storing a collection of 
#' \code{\link{MTuples}} objects. The \code{MTuplesList} class is almost 
#' identical to the \code{\link[GenomicTuples]{GTuplesList}} on which it is 
#' based.
#' 
#' @usage
#' MTuplesList(...)
#' 
#' @details
#'  Please see 
#' \code{\link[GenomicTuples]{GTuplesList}} for a description of available 
#' methods. The only additional methods are \code{\link{methinfo}} and 
#' \code{\link{methtype}}, which are identical to their \code{\link{MTuples}} 
#' counterparts.
#' 
#' @param ... \code{\link{MTuples}} objects. All must contain the same 
#' \code{size} tuples.
#' 
#' @seealso \code{\link{MTuples}}, \code{\link[GenomicTuples]{GTuplesList}}.
#'
#' @export
#' @include MethInfo-class.R
#' @author Peter Hickey
#' @examples 
#' ## TODO
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
    unlistData <- suppressWarnings(do.call("c", unname(listData)))
  }
  
  relist(unlistData, PartitioningByEnd(listData))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @export
setMethod("methinfo", 
          "MTuplesList", 
          function(x) {
            x@unlistData@methinfo
          }
)

#' @export
setMethod("methtype", 
          "MTuplesList", 
          function(x) {
            methtype(x@unlistData@methinfo)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

#' @export
setReplaceMethod("methinfo", 
                 c("MTuplesList", "MethInfo"), 
                 function(x, value) {
                   x@unlistData@methinfo <- value
                   x
                 }
)

#' @export
setReplaceMethod("methtype", 
                 c("MTuplesList", "character"), 
                 function(x, value) {
                   methtype(x@unlistData@methinfo) <- value
                   x
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from MTuples to MTuplesList with extractList() and family.
###
#' @export
setMethod("relistToClass", 
          "MTuples", 
          function(x) {
            "MTuplesList"
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

# Based on GenomicRanges::showList
my_showList <- function(object, showFunction, print.classinfo)
{
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
      showFunction(object[[i]], margin="  ",
                   print.classinfo=print.classinfo)
      if (print.classinfo)
        print.classinfo <- FALSE
      cat("\n")
    }
  } else {
    sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
    if (k >= 3 && cumsumN[3L] <= 20)
      showK <- 3
    else if (k >= 2 && cumsumN[2L] <= 20)
      showK <- 2
    else
      showK <- 1
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
      showFunction(object[[i]], margin="  ",
                   print.classinfo=print.classinfo)
      if (print.classinfo)
        print.classinfo <- FALSE
      cat("\n")
    }
    if (diffK > 0) {
      cat("...\n<", k - showK,
          ifelse(diffK == 1, " more element>\n", " more elements>\n"),
          sep="")
    }
  }
  cat("-------\n")
  cat("seqinfo: ", summary(seqinfo(object)), "\n", sep="")
  cat("methinfo: ", summary(methinfo(object)), "\n", sep = "")
}

#' @export
setMethod("show", 
          "MTuplesList",
          function(object) {
            my_showList(object, showMTuples, FALSE)
          }
)
