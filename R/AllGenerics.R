### -------------------------------------------------------------------------
### methinfo
###

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("methinfo", 
           function(object) {
             standardGeneric("methinfo")
           }
)

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).#' @export
setGeneric("methinfo<-", 
           function(object, value) {
             standardGeneric("methinfo<-")
           }
)

### -------------------------------------------------------------------------
### methtype
###

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("methtype", 
           function(object) {
             standardGeneric("methtype")
           }
)

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("methtype<-", 
           function(object, value) {
             standardGeneric("methtype<-")
           }
)

### -------------------------------------------------------------------------
### methLevel
###

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("methLevel", 
           function(object, ...) {
             standardGeneric("methLevel")
           }
)

### -------------------------------------------------------------------------
### getCoverage
###


# TODO: Decide on a name.
# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("getCoverage", 
           function(object) {
             standardGeneric("getCoverage")
           }
)

### -------------------------------------------------------------------------
### filter
###
###

# TODO: Decide on a better name. Will clash with dplyr's filter()
# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).
#' @export
setGeneric("filter", 
           function(object, ...) {
             standardGeneric("filter")
           }
)
