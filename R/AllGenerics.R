### -------------------------------------------------------------------------
### methinfo
###

#' @export
setGeneric("methinfo", 
           function(object) {
             standardGeneric("methinfo")
           }
)

#' @export
setGeneric("methinfo<-", 
           function(object, value) {
             standardGeneric("methinfo<-")
           }
)

### -------------------------------------------------------------------------
### methtype
###

#' @export
setGeneric("methtype", 
           function(object) {
             standardGeneric("methtype")
           }
)

#' @export
setGeneric("methtype<-", 
           function(object, value) {
             standardGeneric("methtype<-")
           }
)

### -------------------------------------------------------------------------
### methLevel
###

#' @export
setGeneric("methLevel", 
           function(object, ...) {
             standardGeneric("methLevel")
           }
)

### -------------------------------------------------------------------------
### getCoverage
###

## TODO: Decide on a name
#' @export
setGeneric("getCoverage", 
           function(object) {
             standardGeneric("getCoverage")
           }
)
