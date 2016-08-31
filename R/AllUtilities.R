### =========================================================================
### Helper functions not exported.
### For emphasis, __no function in this file will be exported__, because 
### `@keywords internal` applies to the whole file 
### (https://github.com/ramhiser/sparsediscrim/issues/26).
### This means that these functions will not be documented by roxygen2, even
### though the functions have roxygen2 tags.
### =========================================================================

#' Is the methylation type valid
#' 
#' @param methtype A character.
#' @keywords internal
#' @return Returns \code{TRUE} if a valid methylation type, \code{FALSE} 
#' otherwise.
.validMethtype <- function(methtype) {
  if (any(is.na(methtype))) {
    if (length(methtype) == 1L) {
      val <- TRUE
    } else {
      val <- FALSE
    }
  } else {
    VALID_METHTYPE <- c('CG', 'CHG', 'CHH', 'CN')
    val <- all(methtype %in% VALID_METHTYPE)
  }
  
  return(val)
}

#' Create the names of the possible methylation patterns at an m-tuple.
#'
#' This helper function constructs m-tuple names in the correct (alphabetical) 
#' order for a given value of m
#' @param m The size of the m-tuple. Must be an int.
#' 
#' @keywords internal
#' 
#' @examples
#' .makeMethPatNames(1L)
#' .makeMethPatNames(2L)
#' .makeMethPatNames(3L)
#' 
#' @return A character vector
.makeMethPatNames <- function(m) {
  m <- as.integer(m)
  sort(do.call(paste0, 
               expand.grid(lapply(seq_len(m), function(x) {c("M", "U")}))))
}

#' Analogous to rownames and colnames
#' A general version of \code{\link[DSArray]{slicenames}}
#' @keywords internals
.slicenames <- function(x) {
  dimnames(x)[[3L]]
}