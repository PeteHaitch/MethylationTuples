# TODO: See https://stat.ethz.ch/pipermail/bioc-devel/2015-September/008007.html 
# about proper way of writing and documenting S4 generics.

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).

### -------------------------------------------------------------------------
### rowTuples
###

#' @export
setGeneric("rowTuples",
           function(x, ...) {
             standardGeneric("rowTuples")
           }
)

#' @export
setGeneric("rowTuples<-",
           function(x, ..., value) {
             standardGeneric("rowTuples<-")
           }
)


### -------------------------------------------------------------------------
### methinfo
###

#' @export
setGeneric("methinfo", 
           function(x) {
             standardGeneric("methinfo")
           }
)

#' @export
setGeneric("methinfo<-", 
           function(x, value) {
             standardGeneric("methinfo<-")
           }
)

### -------------------------------------------------------------------------
### methtype
###

#' @export
setGeneric("methtype", 
           function(x) {
             standardGeneric("methtype")
           }
)

#' @export
setGeneric("methtype<-", 
           function(x, value) {
             standardGeneric("methtype<-")
           }
)

# TODO: methLevel(), getCoverage() (or similarly named generics)
