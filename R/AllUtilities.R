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
