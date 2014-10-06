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
.make_methpat_names <- function(m){
  if (!is.integer(m) || m < 1){
    stop("'m' must be an int. 'm' must be greater than 0.")
  }
  sort(do.call(paste0, 
               expand.grid(lapply(seq_len(m), function(x) {c('M', 'U')}))))
}