### =========================================================================
### MethPat objects: methylation patterns at genomic tuples
### -------------------------------------------------------------------------
###

## TODO: Make a constructor that is tailored to the MethPat object and not just 
## a thin-wrapper of SummarizedExperiment?
## TODO: Look into using a file-based storage system like NetCDF or the ff 
## package (see ?SummarizedExperiment and 
## https://stat.ethz.ch/pipermail/bioc-devel/2015-September/007992.html).
## TODO: Note, new("MethPat") won't return a valid object although MethPath()
## will. This isn't ideal - I find it merely an annoyance but it may be a 
## bigger problem than I realise.
## TODO: Usage section (will differ from RangedSummarizedExperiment usage 
## section)
#' MethPat instances
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
#' (almost) all the methods documented in 
#' \code{?\link[SummarizedExperiment]{RangedSummarizedExperiment}} also work on 
#' a MethPat object. The methods documented below are additional methods that 
#' are specific to MethPat objects. The key differences are:
#' \itemize{
#'  \item The \code{rowRanges} must be a \code{\link{MTuples}} 
#'  object rather than a \code{\link[GenomicRanges]{GRanges}} object.
#'  \item Certain \code{assays} are required. See \code{assays} argument below.
#' }
#' 
#' @param assays A \code{\link[base]{list}} or 
#' \code{\link[S4Vectors]{SimpleList}} of matrix elements. All elements of the 
#' list must have the same dimensions, and dimension names (if present) must be 
#' consistent across elements and with row names of \code{rowRanges} and 
#' \code{colData}. Specifically, for a MethPat object containing the 
#' methylation patterns at genomic tuples of \code{\link[GenomicTuples]{size}} 
#' \eqn{= m}, there are \eqn{2^m} required assays. For example, for 2-tuples 
#' there are 4 required assays that must be named \code{MM}, \code{MU}, 
#' \code{UM} and \code{UU} (\code{M} = methylated, \code{U} = unmethylated).
#' \strong{TODO:} Should the \code{.makeMethPatNames} function be exported 
#' and referenced here?
#' @param rowTuples A \code{\link{MTuples}} instance describing 
#' the genomic tuple of the methylation loci. Names, if present, become the 
#' row names of the MethPat. The length of the \code{\link{MTuples}} 
#' must equal the number of rows of the matrices in \code{assays}.
#' @param colData An optional, but recommended, 
#' \code{\link[S4Vectors]{DataFrame}} describing the samples. Row names, if 
#' present, become the column names of the MethPat object.
#' @param metadata An optional \code{\link[base]{list}} of arbitrary content 
#' describing the overall experiment.
#' @param ... For \code{MethPat}, S4 methods \code{\link[base]{list}} and 
#' \code{\link[base]{matrix}}, arguments identical to those of the 
#' \code{\link[S4Vectors]{SimpleList}} method. 
#' For \code{rowTuples}, ignored.
#' \strong{TODO}: Check whether 
#' this param can be deleted from docs since the documentation is rather 
#' complex (inherited from RangedSummarizedExperiment).
#' @param x A MethPat object. The \code{rowTuples} setter will also accept a 
#' \link[SummarizedExperiment]{SummarizedExperiment0} object and will first 
#' coerce it to a \link[SummarizedExperiment]{RangedSummarizedExperiment} 
#' before it sets \code{value} on it. \strong{TODO}: Check whether this 
#' actually works.
#' @param value A \link{MTuples} or \link{MTuplesList} object.
#' @param subset An expression which, when evaluated in the context of 
#' \code{rowTuples(x)}, is a logical vector indicating elements of rows to 
#' keep: missing values are taken as false.
#' @param select An expression which, when evaluated in the context of 
#' \code{colData(x)}, is a logical vector indicating elements or rows to 
#' keep: missing values are taken as false.
#' 
#' @section Details:
#' The rows of a MethPat object represent tuples (in genomic coordinates) of 
#' interest. The tuples of interest are described by a \link{MTuples} or a 
#' \link{MTuplesList} object, accessible using the \code{rowTuples} function, 
#' described below. The \link{MTuples} and \link{MTuplesList} classes contain 
#' sequence (e.g., chromosome) name, genomic coordinates, and strand 
#' information, along with methylation-type information (e.g., CG or CHG). Each 
#' tuple can be annotated with additional data; this data might be used to 
#' describe the range (e.g., CpG island status) or to summarized results
#' relevant to the tuple. Rows may or may not have row names; they often will 
#' not.
#' 
#' @section Constructor:
#' MethPat instances are constructed using the \code{MethPat} function with 
#' arguments outlined above.
#' 
#' @section Accessors:
#' In the following code snippets, \code{x} is a MethPat instance.
#' 
#' \describe{
#'  \item{\code{rowTuples(x)}, \code{rowTuples(x) <- value}:}{Get or set the row 
#'  data. \code{value} is a \link{MTuples} or \link{MTuplesList} object. Row 
#'  names of \code{value} must be \code{NULL} or consistent with the existing 
#'  row names of \code{x}}.
#' }
#' 
#' @section MTuples compatibility (rowTuples access):
#' Many \link{MTuples} and \link{MTuplesList} operations are 
#' supported on MethPat objects, using \code{rowTuples}.
#' 
#' Supported operations include: \code{\link{compare}}, 
#' \code{\link{duplicated}}, \code{\link{end}}, \code{\link{end<-}}, 
#' \code{\link{granges}}, \code{\link{match}}, \code{\link{mcols}}, 
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
#' \strong{WARNING:} The preferred getter and setter of tuple information are 
#' \code{rowTuples(x)} and \code{rowTuples(x) <- value}, respectively. 
#' \code{rowRanges(x)}/\code{granges(x)} and \code{rowRanges(x) <- value}/
#' \code{granges(x) <- value} have the same effect since \code{rowTuples} 
#' is merely more informatively named alias to \code{rowRanges} for 
#' \code{MethPat} objects.
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
#'  \item{\code{subset(x, subset, select):}{Create a subset of \code{x} using 
#'  an expression \code{subset} referring to columns of \code{rowTuples(x)} 
#'  (including "seqnames", "start", "end", "width", "strand", and 
#'  \code{names(mcols(x))}) and/or \code{select} referring to column names of 
#'  \code{colData(x)}}}.
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
#'  
#'  \item{\code{combine(x, y, ...)}:}{\code{combine} combines objects with 
#'  different tuples (\code{rowTuples}) and different samples (columns in 
#'  \code{assays}) using an "incomplete" union strategy. Please read 
#'  \code{\link[BiocGenerics]{combine}} for the difference between the union 
#'  and intersection strategies; the current method is "incomplete" because 
#'  it requires that the samples (columns in \code{assays}) are distinct 
#'  across \code{x}, \code{y} and \code{...}. This behaviour may change in 
#'  future versions so that data from the same sample that is stored across 
#'  multiple objects can be safely combined.
#'  }
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
#' @export
#' 
#' @aliases MethPat
#' 
#' @examples
#' ## TODO
setClass("MethPat", 
         contains = "RangedSummarizedExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.MethPat.assays <- function(x) {

  m <- x@rowRanges@size
  if (!is.na(m)) {
    # Check assay names and, if assay names are okay, check counts are all >= 0
    an <- assayNames(x)
    if (is.null(an)) {
      return(paste0("Assay names must include all of: '", 
                    paste0(.makeMethPatNames(m), 
                           collapse = "', '"), "'"))
    } else {
      if (any(is.na(match(.makeMethPatNames(m), an)))) {
        return(paste0("Assay names must include all of: '", 
                      paste0(.makeMethPatNames(m), collapse = "', '"), "'"))
      } else {
        # Note from bsseq: benchmarking shows that min(assay()) < 0 is faster 
        # than any(assay() < 0) if it is false
        if (min(sapply(assays(x)[.makeMethPatNames(m)], 
                       min, na.rm = TRUE), na.rm = TRUE) < 0) {
          return(paste0("All counts of methylation patterns (stored in assays ", 
                 "slot) must be non-negative integers."))
        }
      }
    }
  }
  NULL
}

# TODO: Some sort of validity check on assays, e.g., I require them to be 
#       integers.

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

setValidity2("MethPat", .valid.MethPat)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
MethPat <- function(assays = SimpleList(), 
                    rowTuples = MTuples(), 
                    colData = DataFrame(), 
                    metadata = list()) {
  
  #   if (missing(colData) && 0L != length(assays)) {
  #     nms <- colnames(assays[[1]])
  #     if (is.null(nms) && 0L != ncol(assays[[1]])) {
  #       stop("'MethPat' assay colnames must not be NULL.")
  #     }
  #     colData <- DataFrame(row.names = nms)
  #   }
  
  new("MethPat", SummarizedExperiment(assays = assays, 
                                      rowRanges = rowTuples, 
                                      colData = colData,
                                      metadata = metadata))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

# cbind,MethPat-method and rbind,Methpat-method defined via inheritance to 
# cbind,SummarizedExperiment-method and rbind,SummarizedExperiment-method, 
# respectively.

# I also define a combine,MethPat-method using a (incomplete) union strategy.
# Currently requires that colnames are unique across MethPat objects, i.e., you 
# are combining objects that contain distinct samples.

# WARNING: elementMetadata is dropped
# TODO: A general combine,SummarizedExperiment-method would be great, e.g.,
# combine(x, y, ..., nomatch = NA_integer_), that uses a complete union 
# strategy, i.e., properly combines objects containing potentially duplicate 
# samples and rows.
#' @export
setMethod("combine", 
          c("MethPat", "MethPat"), 
          function(x, y, ...) {
            args <- list(x, y, ...)
            rowTuples <- do.call(c, lapply(args, rowTuples))
            # Remove duplicate tuples
            rowTuples <- unique(rowTuples)
            nt <- length(rowTuples)
            colnames <- unlist(lapply(args, colnames))
            if (anyDuplicated(colnames)) {
              stop("Cannot combine 'MethPat' objects with duplicate colnames.")
            }
            colData <- do.call(rbind, lapply(args, colData))
            an <- lapply(args, assayNames)
            if (any(sapply(an, function(x, y) any(is.na(match(x, y))), 
                           y = an[[1]]))) {
              stop("'MethPat' objects must all contain the same assays.")
            }
            
            # TODO: I suspect that there are faster and more efficient ways to 
            # combine the assays
            # Create assays of the correct dimension (fill with NA_integer_)
            assays <- endoapply(
              assays(args[[1]], withDimnames = FALSE), function(i, nt, colnames) {
                matrix(NA_integer_, nrow = nt, ncol = length(colnames), 
                       dimnames = list(NULL, c(colnames(x), colnames(y))))
              }, nt = nt, colnames = colnames)
            
            # Fill the assays with the values
            for (j in seq_along(args)) {
              ol <- findOverlaps(args[[j]], rowTuples, type = "equal")
              for (i in seq_along(assays)) {
                assays[[i]][subjectHits(ol),
                            match(colnames(args[[j]]), colnames)] <- 
                              assays(args[[j]], withDimnames = FALSE)[[i]]
              }
            }
            assays <- SummarizedExperiment::Assays(assays)
            
            # WARNING: elementMetadata is dropped
            # TOOD: Properly combine elementMetadata slot
            elementMetadata <- x@elementMetadata
            elementMetadata@nrows <- length(rowTuples)

            metadata <- do.call(c, lapply(args, metadata))
            
            BiocGenerics:::replaceSlots(args[[1L]], 
                                        rowRanges = rowTuples,
                                        colData = colData, 
                                        assays = assays,
                                        metadata = metadata,
                                        elementMetadata = elementMetadata)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @export
setMethod("rowTuples", "MethPat", 
          function(x, ...) {
            rowRanges(x)
          }
)

# TODO: Why isn't unique,SummarizedExperiment implemented; at least as 
# unique(rowRanges(x))

#' @export 
setMethod("methinfo", "MethPat", 
          function(x) {
            methinfo(rowRanges(x))
          }
)

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

# Defined via inheritance to split,SummarizedExperiment-method, which in turn 
# calls IRanges::splitAsList. Therefore, IRanges must be listed in Imports. 
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

#' @export
setReplaceMethod("rowTuples", "MethPat",
                 function(x, ..., value) {
                   rowRanges(x) <- value
                   x
                 }
)

#' @export
setReplaceMethod("methinfo", "MethPat", 
                 function(x, value) {
                   methinfo(rowTuples(x)) <- value
                   x
                 }
)

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

#' @export
setMethod("size", "MethPat", 
          function(x) {
            size(rowTuples(x))
          }
)

#' @export
setMethod("tuples", "MethPat", 
          function(x) {
            tuples(rowTuples(x))
          }
)

#' @export
setReplaceMethod("tuples", "MethPat", 
                 function(x, value) {
                   tuples(rowTuples(x)) <- value
                   x
                 }
)

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
# 