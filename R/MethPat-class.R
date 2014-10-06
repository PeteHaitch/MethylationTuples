### =========================================================================
### MTuples objects: methylation patterns at genomic tuples
### -------------------------------------------------------------------------
###

## TODO: Make a constructor that is tailed to the MethPat object and not just a 
## thin-wrapper of SummarizedExperiment.
## TODO: Look into using a file-based storage system like NetCDF or the ff 
## package (see ?SummarizedExperiment).
## TODO: Note, new("MethPat") won't return a valid object although MethPath()
## will. This isn't ideal - I find it merely an annoyance but it may be a 
## bigger problem than I realise.

## TODO: Decide whether to export the class definition; 
## see vignette(topic = 'namespace', package = 'roxygen2').
## TODO: Decide if sample names and methylation type are mandatory arguments 
## and, if so, document these.
## TODO: Usage section (will differ from SummarizedExperiment usage section)
#' MethPat instances
#' 
#' @description
#' The \code{MethPat} class is a matrix-like container where rows represent 
#' genomic tuples of interest and columns represent samples (with sample data 
#' summarized as a \code{\link[S4Vectors]{DataFrame-class}}). A 
#' \code{MethPat} object contains the counts of how many times each methylation 
#' pattern is observed for that genomic tuple in each sample. For example, 
#' there are four possible methylation patterns at 2-tuples: \code{MM}, 
#' \code{MU}, \code{UM} and \code{UU}.
#' 
#' The \code{MethPat} class extends the 
#' \code{\link[GenomicRanges]{SummarizedExperiment}} class. The key 
#' differences are:
#' \itemize{
#'  \item The \code{rowData} must be a \code{\link[GenomicTuples]{GTuples}} 
#'  object rather than a \code{\link[GenomicRanges]{GRanges}} object.
#'  \item Certain \code{assays} are required. See \code{assays} argument below.
#' }
#' 
#' @param assays A \code{\link[base]{list}} or 
#' \code{\link[S4Vectors]{SimpleList}} of matrix elements. All elements of the 
#' list must have the same dimensions, and dimension names (if present) must be 
#' consistent across elements and with row names of \code{rowData} and 
#' \code{colData}. Specifically, for a \code{MethPat} object containing the 
#' methylation patterns at genomic tuples of \code{\link[GenomicTuples]{size}} 
#' \eqn{= m}, there are \eqn{2^m} required assays. For example, for 2-tuples 
#' there are 4 required assays that must be named \code{MM}, \code{MU}, 
#' \code{UM} and \code{UU} (\code{M} = methylated, \code{U} = unmethylated).
#' \strong{TODO:} Should the \code{.make_methpat_names} function be exported 
#' and referenced here?
#' @param rowData A \code{\link[GenomicTuples]{GTuples}} instance describing 
#' the genomic tuples. Row names, if present, become the row names of the 
#' \code{MethPat}. The length of the \code{\link[GenomicTuples]{GTuples}} must 
#' equal the number of rows of the matrices in \code{assays}.
#' @param colData An optional, but recommended, 
#' \code{\link[S4Vectors]{DataFrame}} describing the samples. Row names, if 
#' present, become the column names of the \code{MethPat}.
#' @param exptData An optional \code{\link[S4Vectors]{SimpleList}} of arbitrary 
#' content describing the overall experiment.
#' @param ... For \code{MethPat}, S4 methods \code{\link[base]{list}} and 
#' \code{\link[base]{matrix}}, arguments identical to those of the 
#' \code{\link[S4Vectors]{SimpleList}} method.
#' 
#' For \code{assay}, \code{...} may contain \code{withDimnames}, which is 
#' forwarded to \code{assays}.
#' 
#' For \code{cbind}, \code{rbind}, \code{...} contains \code{MethPat} objects 
#' to be combined.
#' 
#' For other accessors, ignored.
#' @param verbose A \code{logical(1)} indicating whether messages about data 
#' coerction during construction should be printed.
#' @param x, object An instance of \code{MethPat}.
#' @param i,j For \code{assay}, \code{assay<-}, \code{i} is an integer or 
#' numeric scalar; see 'Details' for additional constraints.
#' 
#' For \code{[,MethPat}, \code{[,MethPat<-}, \code{i}, \code{j} are instances 
#' that can act to subset the underlying \code{rowData}, \code{colData}, and 
#' \code{\link[base]{matrix}} elements of \code{assays}.
#' 
#' For \code{[[,MethPat}, \code{[[<-MethPat}, \code{i} is a scalar index (e.g. 
#' \code{character(1)}, or \code{integer(1)}) into a column of \code{colData}.
#' @param subset An expression which, when evaluated in the context of 
#' \code{rowData(x)}, is a logical vector indiciating elements or rows to keep: 
#' missing values are taken as false.
#' @param select An expression which, when evaluated in the context of 
#' \code{colData(x)}, is a logical vector indicating elements or rows to keep: 
#' missing values are taken as false. 
#' @param name A symbol representing the name of a column of \code{colData}.
#' @param withDimnames A \code{logical(1)}, indicating whether dimnames should 
#' be applied to extracted assay elements (this argument is ignored for the 
#' setter \code{assays<-}).
#' @param drop A \code{logical(1)}, ignored by these methods.
#' @param value An instance of a class specified in the S4 method signature or 
#' as outlined in 'Details'.
#' @param deparse.level See \code{\link[base]{cbind}} for a description of this 
#' argument.
#' 
#' @section Constructor:
#' Instances are constructed using the \code{MethPat} function with arguments 
#' outlined aboved.
#' 
#' @section Accessors:
#' In the following code snippets, \code{x} is a \code{MethPat} instance.
#' 
#' \describe{
#'  \item{\code{assays(x)}, \code{assays(x) <- value}:}{Get or set the assays. 
#'  \code{value} is a \code{list} or \code{\link[S4Vectors]{SimpleList}}, each 
#'  element of which is a \code{\link[base]{matrix}} with the same dimensions 
#'  as \code{x}.}
#'  \item{\code{assay(x, i)}, assay(x, i) <- value}:}{A conventient alternative 
#'  (to \code{assays(x)[[i]]})}, \code{assays(x)[[i]] <- value)} to get or set 
#'  the \code{i}th (default first) assay element. \code{value} must be a 
#'  \code{\link[base]{matrix}} of the same dimensions as \code{x}, and with 
#'  dimension names \code{NULL} or consistent with those of \code{x}.
#'  \item{\code{rowData(x), rowData(x) <- value}:}{Get or set the row data. 
#'  \code{value} is a \code{\link[GenomicTuples]{GTuples}} instance. Row names 
#'  of \code{value} must be \code{NULL} or consistent with the existing row 
#'  names of \code{x}.}
#'  \item{\code{colData(x)}, \code{colData(x) <- value}:}{Get or set the column 
#'  data. \code{value} is a \code{\link[S4Vectors]{DataFrame}} instance. Row 
#'  names of \code{value} must be \code{NULL} or consistent with the existing 
#'  columns of \code{x}.}
#'  \item{\code{exptData(x)}, \code{exptData(x) <- value}:}{Get or set the 
#'  experiment data. \code{value} is a \code{\link[base]{list}} or 
#'  \code{\link[S4Vectors]{SimpleList}} instance, with arbitrary content.}
#'  \item{\code{dim(x)}:}{Get the dimensions (tuples x samples) of the 
#'  \code{MethPat} object.}
#'  \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set the 
#'  dimension names. \code{value} is usually a list of length 2, containing 
#'  elements that are either \code{NULL} or vectors of appropriate length for 
#'  the corresponding dimension. \code{value} can be \code{NULL}, which removes 
#'  dimension names. This method implies that \code{rownames}, 
#'  \code{rownames<-}, \code{colnames}, and \code{colnames<-} are all 
#'  available.}
#' 
#' @section GTuples compatibility (rowData access):
#' Many \code{\link[GenomicTuples]{GTuples}} operations are supported on 
#' \code{MetPath} and derived instances, using \code{rowData}. 
#' 
#' \strong{WARNING:} The preferred getter/setter of tuple information is 
#' \code{tuples(x)}/\code{tuples(x) <- value}. In short, the use of 
#' \code{granges(x)}, code{ranges(x)}, \code{ranges(x) <- value}, 
#' \code{start(x)}, \code{start(x) <- value}, \code{end(x)}, 
#' \code{end(x) <- value}, \code{width(x)} and \code{width(x) <- value} is 
#' generally not what is really desired or required when working with 
#' \code{MethPat} objects; see \code{\link[GenomicTuples]{GTuples}} for further 
#' discussion. 
#' 
#' Supported operations include: \code{\link[GenomicTuples]{compare}}, 
#' \code{\link[GenomicTuples]{countOverlaps}}, 
#' \code{\link[GenomicTuples]{distance}}, 
#' \code{\link[GenomicTuples]{distanceToNearest}}, 
#' \code{\link[GenomicTuples]{duplicated}}, 
#' \code{\link[GenomicTuples]{end}} (\strong{not recommended}, see above), 
#' \code{\link[GenomicTuples]{end<-}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{findOverlaps}},
#' \code{\link[GenomicTuples]{follow}}, 
#' \code{\link[GenomicTuples]{granges}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{IPD}},
#' \code{\link[GenomicTuples]{match}},
#' \code{\link[GenomicTuples]{mcols}},
#' \code{\link[GenomicTuples]{mcols<-}},
#' \code{\link[GenomicTuples]{nearest}},
#' \code{\link[GenomicTuples]{order}},
#' \code{\link[GenomicTuples]{overlapsAny}},
#' \code{\link[GenomicTuples]{precede}},
#' \code{\link[GenomicTuples]{ranges}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{ranges<-}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{rank}},
#' \code{\link[GenomicTuples]{relistToClass}},
#' \code{\link[GenomicTuples]{restrict}},
#' \code{\link[GenomicTuples]{seqinfo}},
#' \code{\link[GenomicTuples]{seqinfo<-}},
#' \code{\link[GenomicTuples]{seqnames}},
#' \code{\link[GenomicTuples]{shift}},
#' \code{\link[GenomicTuples]{size}},
#' \code{\link[GenomicTuples]{sort}},
#' \code{\link[GenomicTuples]{split}},
#' \code{\link[GenomicTuples]{start}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{start<-}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{strand}},
#' \code{\link[GenomicTuples]{strand<-}},
#' \code{\link[GenomicTuples]{subsetByOverlaps}},
#' \code{\link[GenomicTuples]{tuples}},
#' \code{\link[GenomicTuples]{tuples<-}},
#' \code{\link[GenomicTuples]{width}} (\strong{not recommended}, see above),
#' \code{\link[GenomicTuples]{width<-}} (\strong{not recommended}, see above).
#' 
#' Not all \code{\link[GenomicTuples]{GTuples}} operations are supported, 
#' because they do not make sense for \code{MethPat} objects (e.g., 
#' \code{length}, \code{name}, \code{as.data.frame}, \code{c}, 
#' \code{splitAsList}), involve non-trivial combination or splitting of rows 
#' (e.g., unique), or have not yet been implemented (\code{window}, 
#' \code{window<-}).
#' 
#' @section Subsetting:
#' \describe{
#'  \item{\code{x[i, j], x[i, j] <- value}:}{Create or replace a subset of 
#'  \code{x}. \code{i}, \code{j} can be \code{numeric}, \code{logical}, 
#'  \code{character}, or \code{missing}. \code{value} must be a 
#'  \code{MethPat} instance with dimensions, dimension names, and assay 
#'  elements consistent with the subset \code{x[i, j]} being replaced.}
#'  \item{\code{subset(x, subset, select)}:}{Create a subset of \code{x} using 
#'  an expression \code{subset} referring to columns of \code{rowData(x)} 
#'  (including \code{seqnames}, \code{start}, \code{end}, \code{width}, 
#'  \code{strand}, and \code{names(mcols(x))}) and / or \code{select} referring 
#'  to column names of \code{colData(x)}.}
#' }
#'  
#' Additional subsetting accessors provide convenient access to \code{colData} 
#' columns
#'  \describe{
#'    \item{\code{x$name, x$name <- value}}{Access or replace column 
#'    \code{name} in \code{x}.}
#'    \item{\code{x[[i, ...]], x[[i, ...]] <- value}}{Access or replace column 
#'    \code{i} in \code{x}.}
#'  }
#'
#' @section Combining:
#' In the code snippets below, \code{...} are \code{MethPat} instances to be 
#' combined.
#' \describe{
#'  \item{\code{cbind(...), rbind(...)}:}{\code{cbind} combines objects with 
#'  identical tuples (\code{rowData}) but different samples (columns in 
#'  \code{assays}). The colnames in \code{colData} must match or an error is 
#'  thrown. Duplicate columns of \code{mcols(rowData(MethPat))} must contain 
#'  the same data.
#'  
#'  \code{rbind} combines objects with different tuples (\code{rowData}) and 
#'  the same subjects in (columns in \code{assays}). Duplicate columns of 
#'  \code{colData} must contain the same data.
#'  
#'  \code{exptData} from all objects are combined into a 
#'  \code{\link[S4Vectors]{SimpleList} with no name checking.}
#'  }
#' }
#' 
#' @author Peter Hickey, building on all the real work of Martin Morgan for the 
#' \code{\link[GenomicRanges]{SummarizedExperiment}} class.
#' 
#' @seealso \code{\link[GenomicRanges]{SummarizedExperiment}}
#' 
#' @examples
#' ## TODO
#' 
setClass('MethPat', 
         contains = "SummarizedExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.valid.MethPat.rowData <- function(object) {
  msg <- NULL
  
  if (!is(object@rowData, "GTuples")) {
    msg <- validMsg(msg, paste0("'rowData' slot of a 'MethPat' object must be ", 
                                "a 'GTuples' object."))
  }
  
  return(msg)
}

.valid.MethPat.assays <- function(object) {
  msg <- NULL
  
  m <- object@rowData@size
  if (!is.na(m)) {
    # Check assay names and, if assay names are okay, check counts are all >= 0
    an <- names(object@assays$field("data"))
    if (is.null(an)) {
      msg <- validMsg(msg, paste0("Assay names must include all of: ", 
                                  paste0(.make_methpat_names(m), 
                                         collapse = ', ')))
    } else {
      if (any(is.na(match(.make_methpat_names(m), an)))) {
        msg <- validMsg(msg, paste0("Assay names must include all of: ", 
                                    paste0(.make_methpat_names(m), 
                                           collapse = ', ')))
      } else {
        # Note from bsseq: benchmarking shows that min(assay()) < 0 is faster 
        # than any(assay() < 0) if it is false
        if (min(sapply(object@assays$field("data")[.make_methpat_names(m)], 
                       min, na.rm = TRUE), na.rm = TRUE) < 0) {
          msg <- validMsg(msg, paste0("All counts of methylation patterns ", 
                                      "(stored in assays slot) must be ", 
                                      "non-negative integers."))
        }
      }
    }
  }
}  
.valid.MethPat <- function(object) {
  
  # First need to check that rowData is an GTuples object.
  # Otherwise some of the .valid.CoMeth.* functions won't work
  msg <- .valid.MethPat.rowData(object)
  if (is.null(msg)){
    
    # Include all other .valid.MethPat.* functions in this vector
    msg <- c(.valid.MethPat.assays(object))
  }
  
  if (is.null(msg)) {
    return(TRUE)
  } else{
    return(msg)
  }
}

setValidity2("MethPat", .valid.MethPat)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
MethPat <- function(assays = SimpleList(), rowData = GTuples(), 
                    colData = DataFrame(), exptData = SimpleList(), ..., 
                    verbose = FALSE) {
  
  if (missing(colData) && 0L != length(assays)) {
    nms <- colnames(assays[[1]])
    if (is.null(nms) && 0L != ncol(assays[[1]])) {
      stop("'MethPat' assay colnames must not be NULL.")
    }
    colData <- DataFrame(row.names = nms)
  }
  
  new("MethPat", SummarizedExperiment(assays = assays, rowData = rowData, 
                                      colData = colData, exptData = exptData, 
                                      ..., verbose = verbose))
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


# TODO: Document combine,MethPat-method
# TODO: A general combine,SummarizedExperiment-method would be great, e.g.,
# combine(x, y, ..., nomatch = NA_integer_), that uses a complete union 
# strategy, i.e., properly combines objects containing potentially duplicate 
# samples and rows.
#' @export
setMethod("combine", 
          c("MethPat", "MethPat"), 
          function(x, y, ...) {
            args <- list(x, y, ...)
            rowData <- unique(do.call(c, lapply(args, function(i) {
              slot(i, "rowData")
            })))
            nr <- length(rowData)
            colnames <- unlist(lapply(args, colnames))
            if (anyDuplicated(colnames)) {
              stop("Cannot combine 'MethPat' objects with duplicate colnames.")
            }
            colData <- do.call(rbind, lapply(args, function(i) {
              slot(i, "colData")
            }))
            an <- lapply(args, function(x) {
              names(x@assays$field('data'))
            })
            if (any(sapply(an, function(x, y) any(is.na(match(x, y))), y = an[[1]]))) {
              stop("'MethPat' objects must all contain the same assays.")
            }
              
            # TODO: I suspect that there are faster and more efficient ways to 
            # combine the assays.
            assays <- endoapply(assays(args[[1]]), function(i, nr, colnames) {
              matrix(NA_integer_, nrow = nr, ncol = length(colnames))
            }, nr = nr, colnames = colnames)
            for (i in seq_along(assays)) {
              for (j in seq_along(args)) {
                ol <- findOverlaps(args[[j]], rowData, type = 'equal')
                assays[[i]][subjectHits(ol), 
                            match(colnames(args[[j]]), colnames)] <- 
                  assays(args[[j]])[[i]]
              }
            }
            assays <- GenomicRanges:::.ShallowSimpleListAssays(data = assays)
            exptData <- do.call(c, lapply(args, exptData))
            initialize(args[[1]], assays = assays, rowData = rowData, 
                       colData = colData, exptData = exptData)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

# Almost all via inheritance to SummarizedExperiment or otherwise implemented 
# in Tuples methods
#' @export
setMethod('granges', 
          'MethPat', 
          function(x, use.mcols = FALSE, ...) {
            granges(rowData(x), use.mcols, ...)
          }
)

# TODO: Why isn't unique,SummarizedExperiment implemented; at least as 
# unique(rowData(x))

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
### Tuples methods
###

#' @export
setMethod("size", 
          "MethPat", 
          function(x) {
            size(rowData(x))
          }
)

#' @export
setMethod("tuples", 
          "MethPat", 
          function(x) {
            tuples(rowData(x))
          }
)

#' @export
setReplaceMethod("tuples", 
                 "MethPat", 
                 function(x, value) {
                   tuples(rowData(x)) <- value
                   x
                 }
)

#' @export
setMethod("IPD", 
          "MethPat", 
          function(x) {
            IPD(rowData(x))
          }
)
