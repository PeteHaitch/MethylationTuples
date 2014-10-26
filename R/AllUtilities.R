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
#' @export
#' 
#' @keywords internal
#' 
#' @examples
#' .make_methpat_names(1L)
#' .make_methpat_names(2L)
#' .make_methpat_names(3L)
#' @return A character vector
.make_methpat_names <- function(m) {
  m <- as.integer(m)
  sort(do.call(paste0, 
               expand.grid(lapply(seq_len(m), function(x) {c('M', 'U')}))))
}

#' Is the methylation type valid
#' 
#' @param methtype A character.
#' @export
#' @keywords internal
#' @return Returns \code{TRUE} if a valid methylation type, \code{FALSE} 
#' otherwise.
.valid_methtype <- function(methtype) {
  if (any(is.na(methtype))) {
    if (length(methtype) == 1L) {
      val <- TRUE
    } else {
      val <- FALSE
    }
  } else {
    VALID_METHTYPE <- c('CG', 'CHG', 'CHH', 'CNN')
    val <- all(methtype %in% VALID_METHTYPE)
  }
  
  return(val)
}

# TODO: Documentation and unit tests
#' Compute correlation coefficient and confidence interval.
#' 
#' Basically, a wrapper around \code{\link[stats]{cor.test}}.
#' 
#' @keywords internal
.my_cor <- function(x, y, method = c("pearson", "kendall", "spearman"), 
                    conf.level = 0.95) {
  method <- match.arg(method)
  z <- try(suppressWarnings(cor.test(x, y, method = method, 
                                     conf.level = conf.level)), silent = TRUE)
  if (is(z, "try-error")) {
    val <- list(cor = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_)
  } else {
    val <- list(cor = z$estimate, CI_lower = z$conf.int[1], 
                CI_upper = z$conf.int[2])
    val[sapply(val, is.null)] <- NA_real_
  }
  return(val)
}
