# TODO: See https://stat.ethz.ch/pipermail/bioc-devel/2015-September/008007.html 
# about proper way of writing and documenting S4 generics.

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).

### -------------------------------------------------------------------------
### rowTuples
###

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("rowTuples",
           function(x, ...) {
             standardGeneric("rowTuples")
           }
)

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("rowTuples<-",
           function(x, ..., value) {
             standardGeneric("rowTuples<-")
           }
)


### -------------------------------------------------------------------------
### methinfo
###

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("methinfo", 
           function(x) {
             standardGeneric("methinfo")
           }
)

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("methinfo<-", 
           function(x, value) {
             standardGeneric("methinfo<-")
           }
)

### -------------------------------------------------------------------------
### methtype
###

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("methtype", 
           function(x) {
             standardGeneric("methtype")
           }
)

#' @importFrom methods setGeneric
#' 
#' @export
setGeneric("methtype<-", 
           function(x, value) {
             standardGeneric("methtype<-")
           }
)

# TODO: methLevel(), getCoverage() (or similarly named generics)
