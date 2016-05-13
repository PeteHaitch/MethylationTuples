### =========================================================================
### MethInfo objects: Basic info about the type of methylation loci in an object
### -------------------------------------------------------------------------
###

#' MethInfo objects
#' 
#' @description 
#' A \code{MethInfo} object is an object that contains basic 
#' information about a set of methylation loci. Currently the only attributes 
#' are the type or context of the methylation loci, but more attributes might 
#' be added in the future as the need arises.
#' 
#' @details
#' Typically \code{MethInfo} objects are not used directly but are part of 
#' higher level objects. Those higher level objects will generally provide a 
#' \code{methinfo} accessor for getting/setting their \code{MethInfo} component.
#' 
#' @section Constructor:
#' \describe{
#' \item{}{
#'  \code{MethInfo(methtype)}: Creates a \code{MethInfo} object.
#'  }
#' }
#' 
#' @section Accessor methods:
#' In the code snippets below, \code{x} is a \code{MethInfo} object.
#'  \describe{
#'    \item{}{
#'      \code{methtype(x)}, \code{methtype(x) <- value}:
#'      Get/set the methylation type of \code{x}. \code{value} must be a 
#'      character vector: \code{"CG"} (\emph{i.e.}, CpG), \code{"CHG"}, 
#'      \code{"CHH"}, \code{"CNN"} or some combination of these, e.g., 
#'      \code{c("CG", "CHG")} (\code{NA_character_} is also allowed, 
#'      but not recommended).
#'    }
#'  }
#' @section Combining MethInfo objects:
#' There is no \code{c} method for \code{MethInfo} objects. Rather, a 
#' \code{merge} method is provided. 
#' 
#' In the code snippet below, \code{x} and \code{y} are \code{MethInfo} objects. 
#'  \describe{
#'    \item{}{
#'      \code{merge(x, y)}:
#'      Merge \code{x} and \code{y} into a single \code{MethInfo} object where 
#'      the methylation type is the union of \code{methtype(x)} and 
#'      \code{methtype(y)}. If the \code{methtype} of any object is missing 
#'      (\code{NA_character_}) then the merged \code{methtype} is also missing 
#'      (\code{NA_character}).
#'    }
#' }
#'  
#' @author Peter Hickey
#' @examples
#' x <- MethInfo("CG")
#' y <- MethInfo(c("CHG", "CG"))
#' y # NB: Pretty-prints methylation type as "CG/CHG"
#' methtype(y) # NB: Returns the methylation type as a character vector
#' merge(x, y)
#'  
#' @aliases MethInfo
#' 
#' @importFrom methods setClass
#'  
#' @export
setClass("MethInfo", 
         representation(
           methtype = "character"),
         prototype(
           methtype = NA_character_)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

#' @rdname MethInfo-class
#' @param x A MethInfo object.
#' @importFrom methods setMethod
#' @export
setMethod("methtype", "MethInfo", 
          function(x) {
            x@methtype
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.MethInfo <- function(x) {
  
  # Include all other .valid.MethInfo.* functions in this vector
  c(.valid.MethInfo.methtype(x))
}

.valid.MethInfo.methtype <- function(x) {
  
  if (!.validMethtype(x@methtype)) {
    return(paste0("Invalid 'methtype'. Must be one or more of 'CG', 'CHG', ", 
                  "'CHH' or 'CNN'"))
  }
  
  NULL
}

#' @importFrom S4Vectors setValidity2
setValidity2("MethInfo", .valid.MethInfo)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

#' @rdname MethInfo-class
#' @param methtype A character vector of methylation types: \code{"CG"} 
#' (\emph{i.e.}, CpG), \code{"CHG"}, \code{"CHH"}, \code{"CNN"} or some 
#' combination of these, e.g., \code{c("CG", "CHG")} (\code{NA_character_} is 
#' also allowed, but not recommended).
#' @importFrom methods new
#' @export
MethInfo <- function(methtype = NA_character_) {
  new("MethInfo", methtype = methtype)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

#' @rdname MethInfo-class
#' @inheritParams methype
#' @param value A character vector of methylation types: \code{"CG"} 
#' (\emph{i.e.}, CpG), \code{"CHG"}, \code{"CHH"}, \code{"CNN"} or some 
#' combination of these, e.g., \code{c("CG", "CHG")} (\code{NA_character_} is 
#' also allowed, but not recommended).
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methtype", c("MethInfo", "character"),
                 function(x, value) { 
                   x@methtype <- value
                   x
                 }
) 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show.
###

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors summary
#' @export
setMethod("summary", "MethInfo", 
          function(object, ...) {
            if (all(is.na(methtype(object)))) {
              ans <- "NA"
            } else {
              ans <- paste0(sort(methtype(object)), collapse = "/")
            }
            paste0(ans, " methylation type")
          }
)

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @export
setMethod("show", "MethInfo", 
          function(object) {
            cat(class(object), " object with ", summary(object), "\n", sep = "")
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###

#' @importFrom methods new
.MethInfo.merge <- function(...) {
  
  args <- unname(list(...))
  # Drop any NULLs
  arg_is_null <- sapply(args, is.null)
  if (any(arg_is_null)) {
    args[arg_is_null] <- NULL
  }
  x <- args[[1L]]
  if (!all(sapply(args, is, class(x)))) {
    stop("all arguments must be ", class(x), " objects (or NULLs)")
  }
  methtype <- unlist(lapply(args, methtype))
  if (any(is.na(methtype))) {
    methtype <- NA_character_
  } else {
    methtype <- sort(unique(methtype))
  }
  new("MethInfo", methtype = methtype) 
}

#' @rdname MethInfo-class
#' @param y,object A MethInfo object.
#' @param ... Additional MethInfo objects.
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors merge
#' @export
setMethod("merge", c("MethInfo", "MethInfo"), 
          function(x, y, ...) {
            .MethInfo.merge(x, y, ...)
          }
)

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors merge
#' @export
setMethod("merge", c("MethInfo", "missing"), 
          function(x, y, ...) { 
            .MethInfo.merge(x, ...)
          }
)

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors merge
#' @export
setMethod("merge", c("MethInfo", "NULL"), 
          function(x, y, ...) {
            .MethInfo.merge(x, ...)
          })

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors merge
#' @export
setMethod("merge", c("NULL", "MethInfo"), 
          function(x, y, ...) {
            .MethInfo.merge(y, ...)
          }
)

#' @rdname MethInfo-class
#' @importFrom methods setMethod
#' @importMethodsFrom S4Vectors merge
#' @export
setMethod("merge", c("missing", "MethInfo"), 
          function(x, y, ...) {
            .MethInfo.merge(y, ...)
          }
)
