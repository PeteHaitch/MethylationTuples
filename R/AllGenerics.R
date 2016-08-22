### -------------------------------------------------------------------------
### methinfo
###

#' Accessing/modifying methylation information
#' 
#' A set of generic functions for getting/setting/modifying the methylation 
#' information stored in an object.
#' 
#' @param x The object from/on which to get/set the methylation information.
#' 
#' @details 
#' The \linkS4class{MethInfo} class plays a central role for the functions 
#' described in this man page. For classes that implement it, \code{methinfo} 
#' should return a \linkS4class{MethInfo} object. Examples of containers that 
#' have a \code{methinfo} gettter and setter: \linkS4class{MTuples} and 
#' \linkS4class{MTuplesList} classes.
#' 
#' @rdname methinfo
#' @importFrom methods setGeneric
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
#' 
#' @importFrom methods setGeneric
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
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("methtype", 
           function(x) {
             standardGeneric("methtype")
           }
)

#' @rdname methinfo
#' @inheritParams methinfo<-
#' 
#' @importFrom methods setGeneric
#' @export
setGeneric("methtype<-", 
           function(x, value) {
             standardGeneric("methtype<-")
           }
)

### -------------------------------------------------------------------------
### strandCollapse
###

#' Collapse an \linkS4class{MTuples}-based object across strand
#' 
#' @param x An \linkS4class{MTuples} object or an object with a slot containing 
#' an \linkS4class{MTuples} object for which a \code{strandCollapse()} method 
#' has been defined.
#' @param ... Further arguments to be passed to or from other methods
#' 
#' @details
#' Collapsing an \linkS4class{MTuples}-based object by strand means to flip 
#' tuples on the negative strand to the positive strand and shifting these by 
#' the appropriate offset. 
#' 
#' Only objects with \code{methtype(x) == "CG"} have easily defined behaviour 
#' when collapsing by strand. This is because 'CG' is a palindromic DNA 
#' sequence (it is its own reverse complement). Thus, if 
#' \code{methtype(x) == "CG"}, the tuples on the negative strand of \code{x} 
#' are flipped to the positive strand and shifted to the left by one position, 
#' e.g., 10 -> 9, 33 -> 32. Tuples on the positive strand are unchanged. 
#' Unstranded tuples are not supported and the method should throw an error if 
#' \code{any(strand(x) == "*" == TRUE)}. All tuples in the returned object will 
#' be on the positive strand.
#' 
#' \code{strandCollapse} is an \emph{endomorphism}, it returns an object of the 
#' same class as \code{x}. However, methods may collapse or aggregate 
#' additional slots in any manner they see fit; for an example, see 
#' \code{?\linkS4class{MethPat}}).
#' 
#' @rdname strandCollapse
#' @importFrom methods setGeneric
#' @export
setGeneric("strandCollapse", 
           function(x, ...) {
             standardGeneric("strandCollapse")
           }
)