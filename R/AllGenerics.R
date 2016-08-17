### -------------------------------------------------------------------------
### methtype
###


#' @rdname MethInfo-class
#' @param x The object from/on which to get/set the methylation information.
#' @export
setGeneric("methtype", 
           function(x) {
             standardGeneric("methtype")
           }
)

#' @rdname MethInfo-class
#' @inheritParams methype
#' @param value A character vector of methylation types: \code{"CG"} 
#' (\emph{i.e.}, CpG), \code{"CHG"}, \code{"CHH"}, \code{"CN"} or some 
#' combination of these, e.g., \code{c("CG", "CHG")} (\code{NA_character_} is 
#' also allowed, but not recommended).
#' @export
setGeneric("methtype<-", 
           function(x, value) {
             standardGeneric("methtype<-")
           }
)
