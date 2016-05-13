### =========================================================================
### Helper functions not exported.
### For emphasis, __no function in this file will be exported__, because 
### `@keywords internal` applies to the whole file 
### (https://github.com/ramhiser/sparsediscrim/issues/26).
### This means that these functions will not be documented by roxygen2, even
### though the functions have roxygen2 tags.
### =========================================================================


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

#' Is the methylation type valid
#' 
#' @param methtype A character.
#' 
#' @keywords internal
#' 
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
    VALID_METHTYPE <- c("CG", "CHG", "CHH", "CNN")
    val <- all(methtype %in% VALID_METHTYPE)
  }
  
  return(val)
}

#' Is an object stranded.
#' 
#' @param x An object with a \code{strand}.
#' 
#' @importFrom S4Vectors runValue
#' 
#' @return \code{TRUE} if object is stranded, \code{FALSE} if unstranded and 
#' \code{NA} if ambiguous.
#' 
#' @keywords internal
.isStranded <- function(x) {
  s <- try(strand(x), silent = TRUE)
  if (is(s, "try-error")) {
    return(FALSE)
  } else {
    levels <- as.vector(unique(runValue(strand(x))))
    if (identical(levels, '*')) {
      return(FALSE)
    } else if ("*" %in% levels && ("+" %in% levels || "-" %in% levels) ||
               length(levels) == 0) {
      return(FALSE)
    } else {
      TRUE
    }
  }
}

#' The slice names of a 3-dimensional array-like object
#' 
#' A non-exported helper function that generalises 
#' \code{\link[DSArray]{slicenames}}.
#' 
#' @param x A 3-dimensional array-like object.
#' 
#' @note Assumes that \code{x} has a well-defined \code{dimnames} getter
#' 
#' @seealso \code{\link[DSArray]{slicenames}}
#' 
#' @importMethodsFrom DSArray slicenames
#'
#' @keywords internal 
.slicenames <- function(x) {
  if (is(x, "DSArray")) {
    return(slicenames(x))
  } else {
    dimnames(x)[[3L]]
  }
}

#' The slice names of a 3-dimensional array-like object
#' 
#' A non-exported helper function that generalises 
#' \code{\link[DSArray]{slicenames<-}}.
#' 
#' @param A 3-dimensional array-like object
#' 
#' @note Assumes that \code{x} has a well-defined \code{dimnames<-} setter
#' 
#' @seealso \code{\link[DSArray]{slicenames<-}}
#' 
#' @importMethodsFrom DSArray slicenames<-
#'
#' @keywords internal 
`.slicenames<-` <- function(x, value) {
  if (is(x, "DSArray")) {
    slicenames(x) <- value
  } else {
    dimnames(x)[[3L]] <- value
  }
  x
}