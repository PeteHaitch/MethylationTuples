### -------------------------------------------------------------------------
### methinfo
###

#' Accessing/modifying methylation information
#' 
#' A set of generic functions for getting/setting/modifying the methylation 
#' information stored in an object.
#' 
#' @details 
#' The \linkS4class{MethInfo} class plays a central role for the functions 
#' described in this man page. For classes that implement it, \code{methinfo} 
#' should return a \linkS4class{MethInfo} object. Examples of containers that 
#' have a \code{methinfo} gettter and setter: \linkS4class{MTuples} and 
#' \linkS4class{MTuplesList} classes.
#' 
#' @rdname methinfo
#' @param x The object from/on which to get/set the methylation information.
#' @export
setGeneric("methinfo", 
           function(x) {
             standardGeneric("methinfo")
           }
)

#' @rdname methinfo
#' @inheritParams methinfo
#' @param value A character vector of methylation types: \code{"CG"} 
#' (\emph{i.e.}, CpG), \code{"CHG"}, \code{"CHH"}, \code{"CN"} or some 
#' combination of these, e.g., \code{c("CG", "CHG")} (\code{NA_character_} is 
#' also allowed, but not recommended).
#' @export
setGeneric("methinfo<-", 
           function(x, value) {
             standardGeneric("methinfo<-")
           }
)

### -------------------------------------------------------------------------
### methtype
###


#' @rdname methinfo
#' @inheritParams methinfo
#' @export
setGeneric("methtype", 
           function(x) {
             standardGeneric("methtype")
           }
)

#' @rdname methinfo
#' @inheritParams methinfo<-
#' @export
setGeneric("methtype<-", 
           function(x, value) {
             standardGeneric("methtype<-")
           }
)


