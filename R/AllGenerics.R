# TODO: See https://stat.ethz.ch/pipermail/bioc-devel/2015-September/008007.html 
# about proper way of writing and documenting S4 generics.

# TODO: Document generic rather than method (consistent with 
# http://r-pkgs.had.co.nz/man.html).

### -------------------------------------------------------------------------
### MethPat
###

#' @rdname MethPat-class
#' @importFrom methods setGeneric
#' @export
setGeneric("MethPat", function(assays, ...) standardGeneric("MethPat"))

### -------------------------------------------------------------------------
### rowTuples
###

#' @rdname MethPat-class
#' @importFrom methods setGeneric
#' @export
setGeneric("rowTuples", function(x, ...) standardGeneric("rowTuples"))

#' @rdname MethPat-class
#' @importFrom methods setGeneric
#' @export
setGeneric("rowTuples<-", function(x, ..., value) standardGeneric("rowTuples<-"))

### -------------------------------------------------------------------------
### methinfo
###

#' Accessing/modifying methylation information
#' 
#' @description A set of generic functions for getting/setting/modifying the 
#' methylation information stored in an object.
#' @param x The object from/on which to get/set the methylation information.
#' @details The \linkS4class{MethInfo} class plays a central role for the 
#' functions described in this man page because:
#' \itemize{
#'  \item All these functions work on a \linkS4class{MethInfo} object.
#'  \item Default \code{methtype} getters and setters are provided. By default, 
#'  \code{methtype(x)} does \code{methtype{methinfo(x)}}. So any class with a 
#'  \code{methinfo} getter will have all the above getters work out-of-the-box. 
#' }
#' 
#' If, in addition, the class defines a \code{methinfo} setter, then all the 
#' corresponding setters also work out-of-the-box.
#' 
#' Examples of containers that have a \code{methinfo} getter and setter: 
#' \linkS4class{MTuples}, \linkS4class{MTuplesList}, and \linkS4class{MethPat} 
#' classes.
#' 
#' @author Peter Hickey
#' @examples 
#' ## TODO
#' @importFrom methods setGeneric
#' @export
setGeneric("methinfo", function(x) standardGeneric("methinfo"))

#' @rdname methinfo
#' @param value Typically a \linkS4class{MethInfo} object for the 
#' \code{methinfo} setter.
#' @importFrom methods setGeneric
#' @export
setGeneric("methinfo<-", function(x, value) standardGeneric("methinfo<-"))

### -------------------------------------------------------------------------
### methtype
###

#' @inheritParams methinfo
#' @rdname methinfo
#' @importFrom methods setGeneric
#' @export
setGeneric("methtype", function(x) standardGeneric("methtype"))

#' @rdname methinfo
#' @importFrom methods setGeneric
#' @inheritParams methinfo<-
#' Either a named or unnamed character vector for the \code{methtype} setter.
#' @export
setGeneric("methtype<-", function(x, value) standardGeneric("methtype<-"))

# TODO: methLevel(), getCoverage() (or similarly named generics)

