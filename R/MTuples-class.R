### =========================================================================
### MTuples objects: genomic tuples of methylation loci
### -------------------------------------------------------------------------
###

#' MTuples objects
#' 
#' @description
#' The \code{MTuples} class is a container for the genomic locations 
#' of tuples of methylation loci in a genome. A CpG is the most common example 
#' of a methylation locus, but non-CG (CHG, CHH, and CN) methylation loci are 
#' also supported.
#' 
#' @details
#' The \code{MTuples} class extends the \linkS4class{GTuples} class. Briefly, 
#' a \linkS4class{GTuples} is itself is an extension of the 
#' \linkS4class{GRanges} class for storing 
#' \emph{genomic tuples} rather than \emph{genomic ranges}. 
#' 
#' The MTuples API will be familiar to anyone who has used the 
#' \linkS4class{GRanges} class from the \pkg{GenomicRanges} package. As shown 
#' in the examples below, an MTuples object behaves much like a 
#' \linkS4class{GRanges} object. The  main difference is in how genomic 
#' overlaps are computed since we are dealing with tuples rather than ranges; 
#' please see \code{?\linkS4class{GTuples}} for full details. Compared to a 
#' \linkS4class{GTuples} object, an MTuples object has an additional slot, 
#' \code{methinfo}, that stores the information about the type 
#' of methylation loci that are stored in the MTuples object. 
#' 
#' @seealso 
#' \itemize{
#'  \item \code{\link{MTuples}}, \code{\link{MTuplesFromGTuples}}, and 
#'    \code{\link{MTuplesFromBSgenome}} for ways to construct an MTuples object
#'  \item \linkS4class{GTuples} for the class from which MTuples inherits
#' }
#' 
#' @examples 
#' #------------------------------------------------------------------------------
#' # Construct two simple MTuples object
#' mt_pos <- MTuples(seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
#'                   tuples = matrix(c(10, 15, 100, 110, 150,
#'                                     15, 27, 110, 150, 154,
#'                                     27, 33, 150, 154, 166), 
#'                                   ncol = 3),
#'                   strand = "+",
#'                   methinfo = MethInfo("CG"))
#' mt_pos
#' 
#' # NOTE: mt_neg stores the cytosine on the negative strand for each CpG in 
#' #       mt_pos
#' mt_neg <- MTuples(seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
#'                   tuples = matrix(c(11, 16, 101, 111, 151,
#'                                     16, 28, 111, 151, 155,
#'                                     28, 34, 151, 155, 167), 
#'                                   ncol = 3),
#'                   strand = "-",
#'                   methinfo = MethInfo("CG"))
#' mt_neg
#' 
#' #------------------------------------------------------------------------------
#' # Getters and setters
#' 
#' # Common getters and setters (inherited from GRanges)
#' seqnames(mt_pos)
#' strand(mt_pos)
#' seqinfo(mt_pos)
#' strand(mt_pos)
#' 
#' # Common getters and setters (inherited from GTuples)
#' tuples(mt_pos)
#' 
#' # New getters and setters defined specfically for MTuples
#' methinfo(mt_pos)
#' methtype(mt_pos)
#' # Change the methinfo
#' methinfo(mt_pos) <- MethInfo(c("CG", "CHG"))
#' methinfo(mt_pos)
#' # And change methtype back
#' methtype(mt_pos) <- "CG"
#' methinfo(mt_pos)
#' 
#' #------------------------------------------------------------------------------
#' # Combining
#' c(mt_pos, mt_neg)
#' 
#' #------------------------------------------------------------------------------
#' # Collapse by strand (currently only works with CG methtype)
#' strandCollapse(mt_neg)
#' 
#' #------------------------------------------------------------------------------
#' # Finding overlaps
#' 
#' # No 'equal' overlaps, even when ignoring strand
#' findOverlaps(mt_pos, mt_neg, type = "equal", ignore.strand = TRUE)
#' # However, there are 'equal' overlaps if we first collapse the negative 
#' # strand tuples onto the foward strand
#' findOverlaps(mt_pos, strandCollapse(mt_neg), type = "equal")
#' 
#' @include MethInfo-class.R
#' @importClassesFrom GenomicTuples GTuples
#' @importFrom methods setClass
#' @export
setClass("MTuples",
         contains = "GTuples",
         representation(
           methinfo = "MethInfo"),
         prototype(
           methinfo = MethInfo())
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

# None as yet. Validity methods "inherited" from GTuples and MethInfo classes.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' MTuples constructors
#' 
#' There are three ways to construct a \linkS4class{MTuples} object:
#' \enumerate{
#'  \item 'From scratch' using \code{MTuples()}
#'  \item From an existing \linkS4class{GTuples} object using 
#'    \code{MTuplesFromGTuples()}
#'  \item By finding and extracting these from a 
#'    \linkS4class{BSgenome} object using \code{MTuplesFromBSgenome()}
#' }
#' 
#' @param seqnames \linkS4class{Rle} object, character vector, or 
#' factor containing the sequence names.
#' @param tuples \link[base]{matrix} object containing the positions of the 
#' tuples. The first column should refer to pos1, the second to pos2, etc.
#' @param strand \linkS4class{Rle} object, character vector, or 
#' factor containing the strand information.
#' @param ... \itemize{
#'  \item For \code{MTuples()}, optional metadata columns. These columns cannot 
#'  be named "\code{start}", "\code{end}", "\code{width}", or "\code{element}".
#'  \item For \code{MTuplesFromBSgenome()}, arguments passed to 
#'  \linkS4class{BSParams}, see 'Extracting MTuples from a BSgenome object'.
#' }  
#' @param seqlengths \code{NULL}, or an integer vector named with 
#' \code{levels(seqnames)} and containing the lengths (or \code{NA}) for each 
#' level in \code{levels(seqnames)}.
#' @param seqinfo \code{NULL}, or a \linkS4class{Seqinfo} object 
#' containing allowed sequence names, lengths (or \code{NA}), and circularity 
#' flag, for each level in \code{levels(seqnames)}. 
#' @param methinfo A \linkS4class{MethInfo} object containing information 
#' about the methylation loci present in the \linkS4class{MTuples} object.
#' 
#' @return An \linkS4class{MTuples} object.
#' 
#' @seealso The \linkS4class{MTuples} class description
#' 
#' @rdname MTuples
#' @importFrom S4Vectors Rle
#' @importFrom GenomicTuples GTuples
#' @importFrom methods new
#' 
#' @examples 
#' # TODO: Examples for MTuples() and MTuplesFromGTuples()
#' 
#' \dontrun{
#' # Find all CG 2-tuples on chr9 of the human reference genome (hg19)
#' # TODO: Should the example use library or require?
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19,
#'                     methinfo = MethInfo("CG"), 
#'                     size = 2,
#'                     exclude = c(paste0("chr", 1:8), "X", "Y", "M", "_"))
#'                     
#' # NOTE: exclude uses regular expressions, so exclude = "chr1" will also 
#'         exclude "chr1", "chr10", "chr11", etc.
#' x <- MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19,
#'                          methinfo = MethInfo("CG"),
#'                          size = 1,
#'                          exclude = "chr1")
#' unique(seqnames(x))
#' }
#' @export
MTuples <- function(seqnames = Rle(), tuples = matrix(), 
                    strand = Rle("*", length(seqnames)), ..., 
                    seqlengths = NULL, seqinfo = NULL, methinfo = MethInfo()) {
  # TODO: Might need to move methinfo ahead of ... to ensure it isn't 
  #       accidentally swallowed up
  gtuples <- GTuples(seqnames = seqnames, tuples = tuples, strand = strand, 
                     ..., seqlengths = seqlengths, seqinfo = seqinfo)
  MTuplesFromGTuples(gtuples, methinfo = methinfo)
}

#' @rdname MTuples
#' @param gtuples A \linkS4class{GTuples} object containing the positions of 
#' the methylation loci as genomic tuples.
#' @inheritParams MTuples
#' @importFrom GenomicTuples GTuples
#' @importFrom methods new
#' 
#' @export
MTuplesFromGTuples <- function(gtuples = GTuples(), methinfo = MethInfo()) {
  if (!is(gtuples, "GTuples")) {
    stop("'gtuples' must be a 'GTuples' object")
  }
  new("MTuples", gtuples, methinfo = methinfo)
}

#' One of the workhorses behind MTuplesFromBSgenome
#' 
#' @param subject A \linkS4class{DNAString} object, usually a sequence 
#' extracted from a \link{BSGenome} object.
#' @inheritParams MTuplesFromBSGenome 
#' @param methtype A character vector of methylation types: \code{"CG"} 
#' (\emph{i.e.}, CpG), \code{"CHG"}, \code{"CHH"} or some 
#' combination of these, e.g., \code{c("CG", "CHG")} (, \code{"CN"} and 
#' \code{NA_character_} are not currently supported). Usually the result of 
#' calling \code{\link{methtype}()} on a \linkS4class{MethInfo} object.
#' 
#' @details This uses \code{fixed = subject} in the call to \code{matchPDict}. 
#' Consequently, for example, \code{methtype = "CHG"} will match the following 
#' sequences in the \code{subject}: the IUPAC-based matches \code{CAG}, 
#' \code{CCG}, and \code{CTG}, as well as the literal \code{CHG}.
#' 
#' @return An \linkS4class{IRanges} object. All ranges have width 1 and point 
#' to the cytosine in the  \code{methtype}, the strand is recorded as 
#' metadata accessible with \code{\link[S4Vectors]{mcols}()}. For example, 
#' a CHH on the reverse strand with foward-strand co-ordinates chr1:3-5 will be 
#' recorded as \code{iranges <- IRanges(5, 5)} with the strand accessed via  
#' \code{mcols(iranges)$strand} and returning the equivalent of 
#' \code{strand(Rle("-"))}.
#' 
#' @importFrom BiocGenerics strand
#' @importFrom Biostrings DNAStringSet matchPDict reverseComplement
#' @importFrom IRanges resize
#' @importFrom S4Vectors DataFrame mcols Rle
.OneTuplesFromDNAString <- function(subject, methtype) {
  
  # Construct the DNAStringSet from methtype
  methtype <- DNAStringSet(methtype)
  n_methtype <- length(methtype)
  # NOTE: According to ?reverseComplement, applying reverseComplement() to the 
  # pattern before calling matchPattern() [matchPDict() in this case] is the 
  # recommended way of searching hits on the minus strand of a chromosome.
  # NOTE: Could be made more efficient by only searching once for palindromic 
  #       methtype (e.g. CG) and making the appropriate transformation of the 
  #       ranges, but this increases the code complexity.
  rc_methtype <- reverseComplement(methtype)
  dnastringset <- c(methtype, reverseComplement(methtype))
  fix <- c(rep("start", n_methtype), 
           rep("end", n_methtype))
  
  mindex <- matchPDict(pdict = dnastringset, 
                       subject = subject,
                       fixed = "subject")
  # Resize mindex and construct the IRanges
  iranges <- do.call(c, mapply(function(mindex_, fix_) {
    resize(mindex_, width = 1L, fix = fix_)
  }, mindex_ = mindex, fix_ = fix))
  # Add strand as metadata
  mcols(iranges) <- DataFrame(strand = strand(Rle(c(rep("+", n_methtype), 
                                                    rep("-", n_methtype)), 
                                                  lengths(mindex))))
  iranges
}

# TODO: Consider exporting this function
#' Another of the workhorses behind MTuplesFromBSgenome
#' 
#' Takes a \linkS4class{MTuples} object containing 1-tuples and constructs all 
#' adjacent m-tuples where m = \code{size}.
#' 
#' @param one_tuples An \linkS4class{MTuples} object containing 1-tuples.
#' @param size An \code{integer(1)} specifying the size of the tuples to be 
#' created.
#' @param ignore.strand When set to \code{TRUE}, the strand information is 
#' ignored in the construction of the m-tuples
#' 
#' @return An \linkS4class{MTuples} object. The output is sorted according to 
#' \code{seqinfo(one_tuples)}.
#' 
#' @importFrom BiocGenerics is.unsorted sort start strand unstrand
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom S4Vectors mcols mcols<- Reduce
.MTuplesFromOneTuples <- function(one_tuples, size, ignore.strand = FALSE) {
  if (size < 1L) {
    stop("'size' must be a positive integer")
  }
  if (size(one_tuples) != 1L || !is(one_tuples, "MTuples")) {
    stop("'one_tuples' must be an 'MTuples' object with 'size(one_tuples)' = 1")
  }
  if (ignore.strand) {
    one_tuples <- unstrand(one_tuples)
  }
  if (size == 1) {
    # Nothing to do except sort
    if (is.unsorted(one_tuples)) {
      one_tuples <- sort(one_tuples)
    }
    return(sort(one_tuples))
  }
  # Loop over seqlevels and then strand to construct m-tuples (m > 1)
  mtuples_by_seqlevel <- 
    # TODO: Parallelisation via BiocParallel
    lapply(split(one_tuples, seqnames(one_tuples)), function(x) {
      # TODO: Parallelisation via BiocParallel (outer loop more important to 
      #       parallelise)
      mtuples_by_strand <- lapply(split(x, strand(x)), function(xx) {
        if (length(xx) == 0L) {
          return(NULL)
        }
        if (is.unsorted(xx)) {
          xx <- sort(xx)
        }
        n <- length(xx)
        if (size > n) {
          stop("Cannot create tuples of size = ", size, " for seqnames = '", 
               unique(seqnames(xx)), "'")
        }
        start <- start(xx)
        tuples <- do.call(cbind, lapply(seq_len(size), function(i) {
          start[seq.int(i, n - size + i)]
        }))
        MTuples(seqnames = seqnames(xx)[seq.int(1, n - size + 1)],
                tuples = tuples,
                strand = strand(xx)[seq.int(1, n - size + 1)],
                seqinfo = seqinfo(xx),
                methinfo = methinfo(xx))
      })
      mtuples_by_strand <-
        mtuples_by_strand[!vapply(mtuples_by_strand, is.null, logical(1L))]
      # TODO: do.call(c, mtuples_by_strand) ought to work
      Reduce(c, mtuples_by_strand)
    })
  mtuples_by_seqlevel <-
    mtuples_by_seqlevel[!vapply(mtuples_by_seqlevel, is.null, logical(1L))]
  # TODO: do.call(c, mtuples_by_seqlevel) ought to work
  Reduce(c, mtuples_by_seqlevel)
}


#' @rdname MTuples
#' 
#' @param bsgenome A \linkS4class{BSgenome} object.
#' @inheritParams MTuples
#' @param size An \code{integer} specifying the size of tuples for which to 
#' search.
#' 
#' @section Extracting MTuples from a BSgenome object:
#' \code{MTuplesFromBSgenome()} will construct an \linkS4class{MTuples} object 
#' by extracting \emph{adjacent} m-tuples from a \linkS4class{BSgenome} object 
#' containing the full genome sequence of a given organism. The \code{'CN'} 
#' \code{methtype} is not currently supported.
#' 
#' The \code{...} are passed through to a call the construction of a 
#' \linkS4class{BSParams} object; this can be useful, for example, to 
#' exclude certain \code{seqlevels} from the \linkS4class{BSgenome} object, 
#' such as unassigned contigs; see 'Examples' and 
#' \code{?\linkS4class{BSParams}} for further detail.
#'
#' @importFrom Biostrings DNAStringSet matchPattern
#' @importFrom BSgenome bsapply
#' @importFrom GenomeInfoDb seqinfo
#' @export
MTuplesFromBSgenome <- function(bsgenome, methinfo, size, ...) {
  # Argument checks
  if (!is(bsgenome, "BSgenome")) {
    stop("'bsgenome' must be a 'BSgenome' object")
  }
  if (!is(methinfo, "MethInfo")) {
    stop("'methinfo' must be a 'MethInfo' object")
  }
  size <- as.integer(size)
  if (size < 1L) {
    stop("'size' must be a positive integer")
  }
  if ("CN" %in% methtype(methinfo)) {
    # NOTE: Need to decide if CN is going to mean C[A,C,G,T] (IUPAC standard) 
    #       or literal CN (Bismark and methtuple standard); see TODO in 
    #       MethInfo-class.R
    stop("'CN' methtype not yet supported")
  }
  
  # Find all 1-tuples
  bsparams <- new("BSParams", X = bsgenome, FUN = .OneTuplesFromDNAString, ...)
  # TODO: Parallelisation via BiocParallel? Would require that bsapply() used 
  #       bplapply()/bpmapply() rather than sapply()
  list_of_one_tuples <- bsapply(bsparams, methtype = methtype(methinfo))
  
  # Construct MTuples object of 1-tuples
  seqnames <- Rle(names(list_of_one_tuples), lengths(list_of_one_tuples))
  tuples <- matrix(do.call(c, unname(lapply(list_of_one_tuples, start))), 
                   ncol = 1)
  strand <- do.call(c, unname(lapply(list_of_one_tuples, function(x) {
    mcols(x)$strand
  })))
  one_tuples <- MTuples(seqnames = seqnames,
                        tuples = tuples,
                        strand = strand,
                        seqinfo = seqinfo(bsgenome),
                        methinfo = methinfo)
  
  # Construct MTuples object of m-tuples (m = size)
  .MTuplesFromOneTuples(one_tuples = one_tuples, 
                        size = size, 
                        ignore.strand = FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

#' @rdname MTuples-class
#' @param object An MTuples object.
#' @importFrom methods callNextMethod setMethod
#' @export
setMethod("show", "MTuples", 
          function(object) {
            callNextMethod() # nocov start
            margin <- "  "
            cat(margin, "methinfo: ", summary(methinfo(object)), "\n", 
                sep = "") # nocov end
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

#' @rdname MTuples-class
#' @param x An MTuples object.
#' @param ... Additional MTuples objects.
#' @param ignore.mcols,recursive See 'Splitting and Combining' section of
#' \code{?}\linkS4class{GTuples}.
#' @importFrom methods setMethod
#' @importFrom GenomicTuples size
#' @export
setMethod("c", "MTuples",
          function(x, ..., ignore.mcols = FALSE, recursive = FALSE) {
            args <- unname(list(x, ...))
            if (!GenomicTuples:::.zero_range(vapply(args, size, integer(1))) && 
                !isTRUE(all(is.na(vapply(args, size, integer(1)))))) {
              stop("Cannot combine ", 
                   paste0(unique(vapply(args, class, character(1))), 
                          collapse = " and "), 
                   " containing tuples of different 'size'.")
            }
            # NOTE: Can't use callNextMethod() because it doesn't support 
            #       ignore.mcols
            val <- GenomicTuples:::.unlist_list_of_GTuples(args,
                                                           ignore.mcols = ignore.mcols)
            methinfo(val) <- do.call(merge, lapply(list(x, ...), methinfo))
            val
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @rdname MTuples-class
#' @inheritParams methinfo
#' @importFrom methods setMethod
#' @export
setMethod("methinfo", "MTuples", 
          function(x) {
            x@methinfo
          }
)

#' @rdname MTuples-class
#' @inheritParams methtype
#' @importFrom methods setMethod
#' @export
setMethod("methtype", "MTuples", 
          function(x) {
            methtype(x@methinfo)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

#' @rdname MTuples-class
#' @inheritParams methinfo<-
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methinfo", c("MTuples", "MethInfo"), 
                 function(x, value) {
                   x@methinfo <- value
                   x
                 }
)

#' @rdname MTuples-class
#' @inheritParams methtype<-
#' @importFrom methods setReplaceMethod
#' @export
setReplaceMethod("methtype", c("MTuples", "character"), 
                 function(x, value) {
                   methtype(x@methinfo) <- value
                   x
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###

#' @importMethodsFrom S4Vectors split
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### strandCollapse
###

#' @rdname MTuples-class
#' @inheritParams strandCollapse
#' @importFrom BiocGenerics strand strand<- unique
#' @importFrom IRanges shift
#' @importFrom methods setMethod
#' @export
setMethod("strandCollapse", "MTuples",
          function(x, ...) {
            if (methtype(x) != "CG") {
              stop("strandCollapse() only supports '", class(x), "' objects ", 
                   "with the 'CG' methtype")
            }
            if (any(strand(x) == "*")) {
              stop("Object contains unstranded tuples")
            }
            x_neg <- x[strand(x) == "-"]
            strand(x_neg) <- "+"
            x_neg <- shift(x_neg, shift = -1L)
            unique(c(x[strand(x) == "+"], x_neg))
          })
