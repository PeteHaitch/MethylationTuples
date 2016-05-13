### =========================================================================
### MethPat objects: methylation patterns at genomic tuples
### -------------------------------------------------------------------------
###

## TODO: Usage section (will differ from RangedSummarizedExperiment usage 
## section)
#' MethPat objects
#' 
#' @description
#' The MethPat class is a matrix-like container where rows represent genomic 
#' tuples of interest (as a \code{MTuples} object) and columns represent 
#' samples (with sample data summarized as a 
#' \code{\link{DataFrame}}). A MethPat object contains the 
#' counts of how many times each methylation pattern is observed for each 
#' genomic tuple in each sample. For example, there are four possible 
#' methylation patterns at 2-tuples: \code{MM}, \code{MU}, \code{UM} and 
#' \code{UU}.
#' 
#' MethPat is a subclass of the 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} class and, as such, 
#' all the slots documented in 
#' \code{?\link[SummarizedExperiment]{RangedSummarizedExperiment}} also exist 
#' for a MethPat object. The key differences are:
#' \itemize{
#'  \item The \code{rowRanges} must be a \code{\link{MTuples}} 
#'  object rather than a \code{\link[GenomicRanges]{GRanges}} object. 
#'  \item The \code{assays} must include an element named 'counts', 
#'  which is a 3-dimensional array-like object storing the counts of 
#'  methylation patterns at each genomic tuple; each row corresponds to a 
#'  genomic tuple, each column to a sample, and each slice to a particular 
#'  methylation pattern. For example, for a MethPat object containing the 
#'  methylation patterns at 10 2-tuples for 5 samples has a 'counts' assay with 
#'  10 rows, 5 columns, and \eqn{2^3 = 8} slices (named \code{MMM}, \code{MMU}, 
#'  \code{MUM}, \code{MUU}, \code{UMM}, \code{UMU}, \code{UUM}, and \code{UUU}, 
#'  where \code{M} = methylated and \code{U} = unmethylated). The counts can be 
#'  stored using any 3-dimensional array-like object. The 
#'  \linkS4class{DSArray} (in-memory) and \linkS4class{HDF5Array} 
#'  (on-disk) representations are particularly useful compared to using the 
#'  \link[base]{array} representation when the 
#'  \code{\link[GenomicTuples]{size}} of the tuples is > 1. 
#' }
#' Similarly, all the methods documented in 
#' \code{?\link[SummarizedExperiment]{RangedSummarizedExperiment}} also work on 
#' a MethPat object. The methods documented below are additional methods that 
#' are specific to MethPat objects.
#' 
#' @section Details:
#' The rows of a MethPat object represent tuples (in genomic coordinates) of 
#' interest. The tuples of interest are described by a \link{MTuples} or a 
#' \link{MTuplesList} object, accessible using the \code{rowTuples} method, 
#' described below. The \link{MTuples} and \link{MTuplesList} classes contain 
#' sequence (e.g., chromosome) name, genomic coordinates, and strand 
#' information, along with methylation-type information (e.g., CG or CHG). Each 
#' tuple can be annotated with additional data; this data might be used to 
#' describe the range (e.g., CpG island status) or to summarized results
#' relevant to the tuple. Rows may or may not have row names; they often will 
#' not.
#' 
#' @section Constructor:
#' MethPat instances are typically constructed by importing data using 
#' \code{\link{read.methtuple}}, but they can also be created using the 
#' \code{MethPat} constructor using arguments similar to the 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} constructor.
#' 
#' @section Accessors:
#' In the following code snippets, \code{x} is a MethPat instance.
#' 
#' \describe{
#'  \item{\code{rowTuples(x)}, \code{rowTuples(x) <- value}:}{Get or set the row 
#'  data. \code{value} is a \link{MTuples} or \link{MTuplesList} object. Row 
#'  names of \code{value} must be \code{NULL} or consistent with the existing 
#'  row names of \code{x}. These are just aliases for \code{rowRanges} and 
#'  \code{rowRanges<-}, respectively.}
#' }
#' 
#' @section MTuples compatibility (rowTuples access):
#' Many \link{MTuples} and \link{MTuplesList} operations are 
#' supported on MethPat objects, using \code{rowTuples}.
#' 
#' Supported operations include: \code{\link{compare}}, 
#' \code{\link{duplicated}}, \code{\link{end}}, \code{\link{end<-}}, 
#' \code{\link[GTuples]{gtuples}}, \code{\link{match}}, \code{\link{mcols}}, 
#' \code{\link{mcols<-}}, \code{\link{methinfo}}, \code{\link{methtype}}, 
#' \code{\link{order}}, \code{\link{ranges}}, \code{\link{ranges<-}}, 
#' \code{\link{seqinfo}}, \code{\link{seqinfo<-}}, \code{\link{seqnames}},
#' \code{\link{sort}}, \code{\link{start}}, \code{\link{start<-}}, 
#' \code{\link{strand}}, \code{\link{strand<-}}, 
#' \code{\link[GenomicTuples]{tuples}}, \code{\link[GenomicTuples]{strand<-}},
#' \code{\link{width}}, \code{\link{width<-}}.
#' 
#' Since the \link{MTuples} class is a subclass of the 
#' \link[GenomicTuples]{GTuples} class, which is in turn a subclass of the 
#' \link[GenomicRanges]{GRanges} class, there are also methods 
#' compatabile with \link[GenomicTuples]{GTuples}/
#' \link[GenomicRanges]{GRanges}. See also 
#' \code{?\link[SummarizedExperiment]{shift}},
#' \code{?\link[SummarizedExperiment]{isDisjoint}}, 
#' \code{?\link[SummarizedExperiment]{coverage}}, 
#' \code{?\link[SummarizedExperiment]{findOverlaps}}, and 
#' \code{?\link[SummarizedExperiment]{nearest}} for more 
#' \emph{GTuples/GRanges compatibility methods}.
#' 
#' Not all \link{MTuples}/\link[GenomicRanges]{GRanges} 
#' operations are supported, because they do not make sense for MethPat objects 
#' (e.g., \code{name}, \code{as.data.frame}, \code{c}, \code{splitAsList}), 
#' involve non-trivial combination or splitting or rows (e.g., \code{disjoin}, 
#' \code{gaps}, \code{reduce}, \code{unique}) or have not yet been implemented 
#' (\code{Ops}, \code{map}, \code{window}, \code{window<-}).
#' 
#' \strong{WARNING:} The use of \code{ranges(x)}, \code{ranges(x) <- value}, 
#' \code{start(x)}, \code{start(x) <- value}, \code{end(x)}, 
#' \code{end(x) <- value}, \code{width(x)} and \code{width(x) <- value} are 
#' discouraged (although not forbidden) since these are generally not what is 
#' really desired or required when working with MethPat objects; see 
#' \link[GenomicTuples]{GTuples} for further discussion. 
#' 
#' @section Subsetting:
#' In the code snippets below, \code{x} is a MethPat object.
#' 
#' \describe{
#'  \item{\code{subset(x, subset, select):}}{Create a subset of \code{x} using 
#'  an expression \code{subset} referring to columns of \code{rowTuples(x)} 
#'  (including "seqnames", "start", "end", "width", "strand", and 
#'  \code{names(mcols(x))}) and/or \code{select} referring to column names of 
#'  \code{colData(x)}.}
#' }
#' 
#' @section Extension:
#' MethPat is implemented as an S4 class, and can be extended in the usual way, 
#' using \code{contains = "Methpat"} in the new class definition.
#' 
#' @section Combining:
#' In the code snippets below, \code{x}, \code{y} and \code{...} are 
#' MethPat instances to be combined. All \code{MethPat} instances must have the 
#' same \code{\link{size}} tuples and have compatible \link{Seqinfo} and 
#' \link{MethInfo}.
#' 
#' \describe{
#'  \item{\code{cbind(...), rbind(...)}:}{\code{cbind} combines objects with 
#'  identical tuples (\code{rowTuples}) but different samples (columns in 
#'  \code{assays}). The colnames in \code{colData} must match or an error is 
#'  thrown. Duplicate columns of \code{mcols(rowRanges(MethPat))} must 
#'  contain the same data.
#'  
#'  \code{rbind} combines objects with different tuples (\code{rowTuples}) and 
#'  the same subjects in (columns in \code{assays}). Duplicate columns of 
#'  \code{colData} must contain the same data.
#'  
#'  \code{metadata} from all objects are combined into a 
#'  \code{list} with no name checking.}
#' }
#' 
#' @author Peter Hickey, \email{peter.hickey@@gmail.com}, building on all the 
#' real work of Martin Morgan for the 
#' \link[SummarizedExperiment]{RangedSummarizedExperiment} class.
#' 
#' @seealso 
#' \itemize{
#'  \item \link{RangedSummarizedExperiment} objects.
#'  \item \link{SummarizedExperiment0} objects.
#'  \item \link[SummarizedExperiment]{shift}, 
#'  \link[SummarizedExperiment]{isDisjoint}, 
#'  \link[SummarizedExperiment]{coverage}, 
#'  \link[SummarizedExperiment]{findOverlaps}, and 
#'  \link[SummarizedExperiment]{nearest} for more 
#'  \emph{GTuples/GRanges compatibility methods}.
#'  \item \link{MTuples} objects.
#'  \item \link[GenomicTuples]{GTuples} objects in the \pkg{GenomicTuples} 
#'  package.
#' }
#' 
#' @aliases MethPat
#' 
#' @examples
#' ## TODO
#' 
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom methods setClass
#' @export
setClass("MethPat", 
         contains = "RangedSummarizedExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

#' @importFrom methods is
.valid.MethPat.rowTuples <- function(x) {
  if (!is(x@rowRanges, "MTuples")) {
    return(paste0("'rowTuples' of '", class(x), "' object must be an ", 
                  "'MTuples' object."))
  }
  NULL
}


#' @importMethodsFrom DSArray slicenames
#' @importMethodsFrom SummarizedExperiment assayNames assay
.valid.MethPat.assays <- function(x) {

  m <- size(x@rowRanges)
  if (!is.na(m)) {
    # Check assay names and, if assay names are okay, check counts are all >= 0
    an <- assayNames(x)
    if (is.null(an) || !("counts" %in% an)) {
      return("Assays must include an element named 'counts'")
    } else {
      if (any(is.na(match(.makeMethPatNames(m), 
                          .slicenames(assay(x, "counts")))))) {
        return(paste0("'counts' slicenames must be: '", 
                      paste0(.makeMethPatNames(m), collapse = "', '"), "'"))
      } else {
        # Note from bsseq: benchmarking shows that min(assay()) < 0 is faster 
        # than any(assay() < 0) if it is false
        if (min(assay(x, "counts")) < 0) {
          return("All 'counts' must be non-negative integers.")
        }
      }
    }
  }
  NULL
}

# TODO: Some sort of validity check on assays? E.g., they should be counts 
#       (integers) but this could be a costly check when the assays use the 
#       DSArray or HDF5Array class.
.valid.MethPat <- function(x) {
  
  # First need to check that rowTuples is an MTuples object.
  # Otherwise some of the .valid.MethPat.* functions won't work
  msg <- .valid.MethPat.rowTuples(x)
  if (is.null(msg)) {
    
    # Include all other .valid.MethPat.* functions in this vector
    c(.valid.MethPat.assays(x))
  } else {
    msg
  }
}

#' @importFrom S4Vectors setValidity2
setValidity2("MethPat", .valid.MethPat)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @rdname MethPat-class
#' @param assays A \link{list} of \linkS4class{SimpleList} of 
#' assays. Must include a 3-dimensional array-like assay named 'counts'; see 
#' 'Description' for further details.
#' @param rowTuples A \linkS4class{MTuples} object; see 'Description' for 
#' further details.
#' @param ... Additional arguments passed down to the 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} constructor, such 
#' as \code{colData} and \code{metadata}. Note that this may not include 
#' \code{rowRanges} or \code{rowData} since MethPat objects instead use 
#' \code{rowTuples}.
#' @importFrom methods callNextMethod setMethod
#' @importFrom S4Vectors SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment
#' 
#' @export
setMethod("MethPat", "ANY", 
          function(assays, rowTuples, ...) {
            dots <- list(...)
            if ("rowRanges" %in% names(dots) | 
                "rowData" %in% names(dots)) {
              stop("'...' must not include 'rowRanges' or 'rowData'")
            }
            if (!is.list(assays) | !is(assays, "List")) {
              assays <- SimpleList(counts = assays)
            }
            if (missing(rowTuples)) {
              rowTuples <- MTuples()
            }
           new("MethPat", SummarizedExperiment(assays = assays, 
                                               rowRanges = rowTuples, 
                                               ...))
          }
)

#' @rdname MethPat-class
#' @inheritParams MethPat
#' @importFrom methods callNextMethod setMethod
#' @importFrom S4Vectors SimpleList
#' @export
setMethod("MethPat", "missing", 
          function(assays, rowTuples, ...) {
            dots <- list(...)
            if ("rowRanges" %in% names(dots) | 
                "rowData" %in% names(dots)) {
              stop("'...' must not include 'rowRanges' or 'rowData'")
            }
            if (missing(rowTuples)) {
              rowTuples <- MTuples()
            }
            new("MethPat", SummarizedExperiment(SimpleList(), 
                                                rowRanges = rowTuples, 
                                                ...))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

# cbind,MethPat-method and rbind,Methpat-method defined via inheritance to 
# cbind,SummarizedExperiment-method and rbind,SummarizedExperiment-method, 
# respectively.

# TODO: A general combine,SummarizedExperiment-method would be great, e.g.,
#       combine(x, y, ..., nomatch = NA_integer_), that uses a complete union 
#       strategy, i.e., properly combines objects containing potentially 
#       duplicate samples and rows. This would require that any array-like 
#       class that are used as elements in the assays slot have a well-defined 
#       combine() method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @rdname MethPat-class
#' @param x A MethPat object.
#' @importFrom methods setMethod
#' @importMethodsFrom SummarizedExperiment rowRanges
#' @export
setMethod("rowTuples", "MethPat", 
          function(x, ...) {
            rowRanges(x)
          }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @importFrom methods setMethod
#' @importMethodsFrom SummarizedExperiment rowRanges
#' @export 
setMethod("methinfo", "MethPat", 
          function(x) {
            methinfo(rowRanges(x))
          }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @importFrom methods setMethod
#' @importMethodsFrom SummarizedExperiment rowRanges
#' @export
setMethod("methtype", "MethPat", 
          function(x) {
            methtype(rowRanges(x))
          }
)

# TODO: methLevel,MethPat-method, getCoverage,MethPat-method (depending on what
# I call the generic in R/AllGenerics.R)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###

# TODO: split,MethPat,Rle-method has an ambiguous signature. According to 
#       http://adv-r.had.co.nz/S4.html this isn't desirable. This is the note - 
#       "Note: method with signature ‘SummarizedExperiment#ANY’ chosen for 
#       function ‘split’, target signature ‘MethPat#Rle’ "ANY#Vector" would 
#       also be valid."
#       The same warning can be triggered by:
#       example("RangedSummarizedExperiment)
#       split(rse, seqnames(rse))
#       ?RangedSummarizedExperiment does note that split() isn't supported for 
#       RangedSummarizedExperiment objects, but it "appears" to work, which 
#       might be confusing.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

# Most defined via inheritance to SummarizedExperiment or implemented in Tuples 
# methods

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @param value A replacement object of the appropriate class.
#' @importFrom methods setReplaceMethod
#' @importMethodsFrom SummarizedExperiment "rowRanges<-"
#' @export
setReplaceMethod("rowTuples", "MethPat",
                 function(x, ..., value) {
                   rowRanges(x) <- value
                   x
                 }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples<-
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methinfo", "MethPat", 
                 function(x, value) {
                   methinfo(rowTuples(x)) <- value
                   x
                 }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples<-
#' @importFrom methods setReplaceMethod
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methtype", "MethPat", 
                 function(x, value) {
                   methtype(rowTuples(x)) <- value
                   x
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @importFrom methods setMethod
#' @importMethodsFrom GenomicTuples size
#' @export
setMethod("size", "MethPat", 
          function(x) {
            size(rowTuples(x))
          }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @importFrom methods setMethod
#' @importMethodsFrom GenomicTuples tuples
#' @export
setMethod("tuples", "MethPat", 
          function(x) {
            tuples(rowTuples(x))
          }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples<-
#' @importFrom methods setReplaceMethod
#' @importMethodsFrom GenomicTuples "tuples<-"
#' @export
setReplaceMethod("tuples", "MethPat", 
                 function(x, value) {
                   tuples(rowTuples(x)) <- value
                   x
                 }
)

#' @rdname MethPat-class
#' @inheritParams rowTuples
#' @importFrom methods setMethod
#' @importMethodsFrom GenomicTuples IPD
#' @export
setMethod("IPD", "MethPat", 
          function(x) {
            IPD(rowTuples(x))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method
###

# TODO: Re-write based on show,SummarizedExperiment-method. Be sure to include 
# methinfo and rename rowRanges rowTuples
