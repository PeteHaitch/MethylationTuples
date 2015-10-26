### =========================================================================
### SparseMethPat objects: sparse methylation patterns at genomic tuples
### -------------------------------------------------------------------------
###

#' SparseMethPat objects
#'
#' @rdname SparseMethPat
#' 
#' @importFrom methods setClass
#'
#' @export
setClass("SparseMethPat", 
         contains = "RangedSparseSummarizedExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#' @importMethodsFrom SparseSummarizedExperiment sparseAssayNames
.valid.SparseMethPat <- function(x) {
  
  # Check that rowTuples slot is an MTuples object
  .valid.MP.rowTuples(x)
  
  # Check sparseAssays
  # NOTE: Only check these if object has size.
  m <- x@rowRanges@size
  if (!is.na(m)) {
    # Check that there is a sparse assay called 'counts'
    if (!"counts" %in% sparseAssayNames(x)) {
      return("'sparseAssays' slot must include a sparse assay named 'counts'")
    }
    
    # Check colnames of data elements in counts slot
    data_colnames <- lapply(x@sparseAssays[["counts"]], function(sample) {
      colnames(sample[["data"]])
    })
    cn <- .makeMethPatNames(m)
    if (!all(vapply(data_colnames, identical, logical(1), cn))) {
      return(paste0("colnames of 'data' element in 'counts' sparse assay must ", 
                    "be: ", paste0(cn, collapse = ", ")))
    }
    
    # Check all data elements are >= 0. 
    # NOTE: min(y) < 0 is faster than any(y < 0) if FALSE.
    data_min <- lapply(x@sparseAssays[["counts"]], function(sample) {
      min(sample[["data"]])
    })
    data_min <- unlist(data_min, use.names = FALSE)
    if (isTRUE(any(data_min < 0)) || isTRUE(any(is.na(data_min)))) {
      return(paste0("'data' elements in 'counts' sparse assay must be ", 
                    "non-negative and non-NA"))
    }
  }
  
  NULL
}

#' @importFrom S4Vectors setValidity2
setValidity2("SparseMethPat", .valid.SparseMethPat)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

# TODO: Unify the MethPat() and SparseMethPat() constructors, e.g,
# MethPat <- function(mtuples, counts, colData, metadata, sparse), where 
# sparse = TRUE returns a SparseMethPat object, otherwise a MethPat object.

#' @rdname SparseMethPat
#'
#' @importFrom methods new
#' @importFrom SparseSummarizedExperiment sparseAssays 
#'                                        SparseSummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @export
SparseMethPat <- function(sparseAssays = sparseAssays(), 
                          rowTuples = MTuples(), 
                          colData = DataFrame(), 
                          metadata = list()) {

  new("SparseMethPat", SparseSummarizedExperiment(sparseAssays = sparseAssays, 
                                                  rowRanges = rowTuples, 
                                                  colData = colData,
                                                  metadata = metadata))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

# cbind,SparseMethPat-method and rbind,SparseMethpat-method defined via 
# inheritance to 
# cbind,RangedSparseSummarizedExperiment-method and 
# rbind,RangedSparseSummarizedExperiment-method, respectively.

# I also define a combine,MethPat-method using a (incomplete) union strategy.
# Currently requires that colnames are unique across MethPat objects, i.e., you 
# are combining objects that contain distinct samples.

# TODO: A general combine,RangedSparseSummarizedExperiment-method would be 
# great, e.g., combine(x, y, ..., nomatch = NA_integer_):
# (1) [easier] using an incomplete union strategy, combining objects with 
# distinct samples but possibly overlapping tuples.
# (2) [hard] using a complete union strategy, i.e., properly combines objects 
# containing potentially duplicate samples and tuples.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @rdname SparseMethPat
#' 
#' @importFrom methods setMethod
#' @importMethodsFrom SummarizedExperiment rowRanges
#'
#' @export
setMethod("rowTuples", "SparseMethPat", 
          function(x, ...) {
            rowRanges(x)
          }
)

# TODO: Why isn't unique,SummarizedExperiment implemented; at least as 
# unique(rowRanges(x))

#' @rdname SparseMethPat
#'
#' @export
setMethod("methinfo", "SparseMethPat", 
          function(x) {
            methinfo(rowRanges(x))
          }
)

#' @rdname SparseMethPat
#'
#' @importFrom methods setMethod
#' @importMethodsFrom SummarizedExperiment rowRanges
#' 
#' @export
setMethod("methtype", "SparseMethPat", 
          function(x) {
            methtype(rowRanges(x))
          }
)

# TODO: methLevel,SparseMethPat-method, getCoverage,SparseMethPat-method 
# (depending on what I call the generic in R/AllGenerics.R)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###

# Defined via inheritance to split,RangedSummarizedExperiment-method, which in 
# turn calls IRanges::splitAsList. Therefore, IRanges must be listed in Imports. 
# NB: Can't selectivly import IRanges::splitAsList because this function calls 
# other functions listed in IRanges, hence it is easiest to simply import the 
# entire IRanges package.

# TODO: split,MethPat,Rle-method has an ambiguous signature. According to 
# http://adv-r.had.co.nz/S4.html this isn't desirable. This is the note - 
# "Note: method with signature ‘SummarizedExperiment#ANY’ chosen for function 
# ‘split’, target signature ‘MethPat#Rle’ "ANY#Vector" would also be valid."

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

# Most defined via inheritance to SummarizedExperiment or implemented in Tuples 
# methods

#' @rdname SparseMethPat
#' 
#' @importFrom methods setReplaceMethod
#' @importMethodsFrom SummarizedExperiment "rowRanges<-"
#'
#' @export
setReplaceMethod("rowTuples", "SparseMethPat",
                 function(x, ..., value) {
                   rowRanges(x) <- value
                   x
                 }
)

#' @rdname SparseMethPat
#' 
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("methinfo", "SparseMethPat", 
                 function(x, value) {
                   methinfo(rowTuples(x)) <- value
                   x
                 }
)

#' @rdname SparseMethPat
#' 
#' @importFrom methods setReplaceMethod
#'
#' @export
setReplaceMethod("methtype", "SparseMethPat", 
                 function(x, value) {
                   methtype(rowTuples(x)) <- value
                   x
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

#' @rdname SparseMethPat
#' 
#' @importFrom methods setMethod
#' @importMethodsFrom GenomicTuples size
#'
#' @export
setMethod("size", "SparseMethPat", 
          function(x) {
            size(rowTuples(x))
          }
)

#' @rdname SparseMethPat
#' 
#' @importFrom methods setMethod
#' @importMethodsFrom GenomicTuples tuples
#'
#' @export
setMethod("tuples", "SparseMethPat", 
          function(x) {
            tuples(rowTuples(x))
          }
)

#' @rdname SparseMethPat
#' 
#' @importFrom methods setReplaceMethod
#' @importMethodsFrom GenomicTuples "tuples<-"
#'
#' @export
setReplaceMethod("tuples", "SparseMethPat", 
                 function(x, value) {
                   tuples(rowTuples(x)) <- value
                   x
                 }
)

#' @rdname SparseMethPat
#'
#' @importMethodsFrom GenomicTuples IPD
#'
#' @export
setMethod("IPD", "SparseMethPat", 
          function(x) {
            IPD(rowTuples(x))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method
###

# TODO: Re-write based on show,SummarizedExperiment-method. Be sure to include 
# methinfo and rename rowRanges rowTuples
# 
