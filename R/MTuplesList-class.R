### =========================================================================
### MTuplesList objects
### -------------------------------------------------------------------------
###

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
#' Please see \linkS4class{GTuplesList} for a description of available 
#' methods. The only additional methods are \code{methinfo} and 
#' \code{methtype}, which are identical to their \linkS4class{MTuples}
#' counterparts.
#' 
#' @param ... \code{\link{MTuples}} objects. All must contain the same 
#' \code{size} tuples.
#' 
#' @seealso \linkS4class{MTuples}, \linkS4class{GTuplesList}
#' 
#' @aliases MTuplesList
#'
#' @include MethInfo-class.R
#' @examples 
#' ## TODO
#' 
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

#' @importFrom IRanges PartitioningByEnd
#' @importFrom BiocGenerics relist
#' 
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
      stop("all elements in '...' must be 'MTuples' objects")
    }
    if (!GenomicTuples:::.zero_range(sapply(listData, size)) && 
        !isTRUE(all(is.na(sapply(listData, size))))) {
      stop("all 'MTuples' objects in '...' must have the same 'size'")
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
#' @inheritParams methinfo
#' @importFrom methods setMethod
#' 
#' @export
setMethod("methinfo", "MTuplesList", 
          function(x) {
            x@unlistData@methinfo
          }
)

#' @rdname MTuplesList-class
#' @inheritParams methtype
#' @importFrom methods setMethod
#' 
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
#' @inheritParams methinfo<-
#' @importFrom methods setReplaceMethod
#' 
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
#' 
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
#' @importFrom methods setMethod
#' @importMethodsFrom IRanges relistToClass
#' 
#' @export
setMethod("relistToClass", "MTuples", 
          function(x) {
            "MTuplesList"
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

#' @rdname MTuplesList-class
#' @param object An MTuplesList object.
#' @importFrom methods callNextMethod setMethod
#' 
#' @export
setMethod("show", "MTuplesList",
          function(object) {
            callNextMethod() # nocov start
            margin <- ""
            cat(margin, "methinfo: ", summary(methinfo(object)), "\n", 
                sep = "") # nocov end
          }
)
