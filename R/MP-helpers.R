### =========================================================================
### MP helper functions
### -------------------------------------------------------------------------
###
### These are helper functions for methods defined for MethPat and 
### SparseMethPat objects, collectively abbrevaited as MP objects.
###
### None of these functions are exported but their functionality is made
### available via the appropriate S4 method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#' @importFrom methods is
#' @keywords internal
.valid.MP.rowTuples <- function(x) {
  if (!is(x@rowRanges, "MTuples")) {
    return(paste0("'rowTuples' of '", class(x), "' object must be an ", 
                  "'MTuples' object."))
  }
  NULL
}
