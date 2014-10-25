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

### -------------------------------------------------------------------------
### betaVal
###

#' @export
setGeneric("betaVal", 
           function(x, ...) {
             standardGeneric("betaVal")
           }
)

### -------------------------------------------------------------------------
### getCov 
###

## TODO: Decide on a name and export
#' @export
setGeneric("getCoverage", 
           function(x) {
             standardGeneric("getCoverage")
           }
)
