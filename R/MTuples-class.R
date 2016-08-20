### =========================================================================
### MTuples objects: genomic tuples of methylation loci
### -------------------------------------------------------------------------
###

#' MTuples objects
#' 
#' @description
#' The \code{MTuples} class is a container for the genomic locations 
#' of methylation loci in a genome.
#' 
#' @details
#' The \code{MTuples} class extends the \linkS4class{GTuples} class by adding 
#' the \code{methinfo} slot, which records the information about the type of 
#' methylation loci that are stored in the MTuplesList object.
#'
# TODO: Include this paragraph once added
# @seealso \code{\link{findMTuples}} to find genomic tuples of methylation
# loci in a \linkS4class{BSgenome} object, and \linkS4class{GTuples} for the
# class from which MTuples inherits.
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

# None as yet. Validity methods "inherited" from GRanges and MethInfo classes.

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
#' @seealso \linkS4class{MTuples}
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


#' Another of the workhorses behind MTuplesFromBSgenome
#' 
#' @param seqnames A \code{character(1)} used for the \code{seqnames} of the 
#' returned 
#' object.
#' @param one_tuples One element of the output of 
#' \code{.OneTuplesFromDNAString()}.
#' @param size An \code{integer(1)} specifying the size of the tuples to be 
#' created.
#' @param seqinfo A \linkS4class{Seqinfo} object.
#' @param methinfo A \linkS4class{MethInfo} object.
#' @param ignore.strand When set to \code{TRUE}, the strand information is 
#' ignored in the construction of the m-tuples
#' 
#' @return An \linkS4class{MTuples} object.
#' 
#' @importFrom BiocGenerics is.unsorted sort start strand unstrand
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols mcols<-
.MTuplesFromOneTuples <- function(seqnames, one_tuples, size, seqinfo, methinfo,
                                  ignore.strand = FALSE) {
  strand <- mcols(one_tuples)$strand
  mcols(one_tuples) <- NULL
  granges <- GRanges(seqnames = seqnames,
                     ranges = one_tuples,
                     strand = strand)
  if (ignore.strand) {
    granges <- unstrand(granges)
  }
  mtuples_by_strand <- lapply(split(granges, strand(granges)), function(gr) {
    if (length(gr) == 0L) {
      return(NULL)
    }
    if (is.unsorted(gr)) {
      gr <- sort(gr)
    }
    n <- length(gr)
    start <- start(gr)
    tuples <- do.call(cbind, lapply(seq_len(size), function(i) {
      if (size > n) {
        warning("Cannot create tuples of size = ", size)
        return(NULL)
      }
      start[seq.int(i, n - size + i)]
    }))
    MTuples(seqnames = seqnames,
            tuples = tuples,
            strand = strand(gr)[seq.int(1, n - size + 1)],
            seqinfo = seqinfo,
            methinfo = methinfo)
  })
  mtuples_by_strand <- mtuples_by_strand[!vapply(mtuples_by_strand, is.null,
                                                 logical(1L))]
  # TODO: do.call(c, mtuples_by_strand) ought to work
  Reduce(c, mtuples_by_strand)
}


# UP TO HERE: MTuplesFromBSGenome() (replacing findMTuples())
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
#' containing the full genome sequence of a given organism.
#' The \code{...} are passed through to a call the construction of a 
#' \linkS4class{BSParams} object. This can be useful, for example, to 
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
    stop("'bsgenome' must be a 'BSgenome' object.")
  }
  if (!is(methinfo, "MethInfo")) {
    stop("'methinfo' must be a 'MethInfo' object.")
  }
  size <- as.integer(size)
  if (size < 1L) {
    stop("'size' must be a positive integer.")
  }
  if ("CN" %in% methtype(methinfo)) {
    # TODO: Need to decide if CN is going to mean C[A,C,G,T] (IUPAC standard) 
    #       or literal CN (Bismark and methtuple standard); see TODO in 
    #       MethInfo-class.R
    stop("'CN' methtype not yet supported")
  }
  
  # Find all 1-tuples
  bsparams <- new("BSParams", X = bsgenome, FUN = .OneTuplesFromDNAString, ...)
  # TODO: Parallelisation via BiocParallel? Would require that bsapply() used 
  #       bplapply()/bpmapply() rather than sapply()
  list_of_one_tuples <- bsapply(bsparams, methtype = methtype(methinfo))
  
  # UP TO HERE: Construct all m-tuples from 1-tuples
  # TODO: Parallelisation via BiocParallel
  list_of_mtuples <- mapply(.MTuplesFromOneTuples, 
                            seqnames = names(list_of_one_tuples),
                            one_tuples = list_of_one_tuples,
                            MoreArgs = list(size = size,
                                            seqinfo = seqinfo(bsgenome),
                                            methinfo = methinfo,
                                            ignore.strand = FALSE))
  list_of_mtuples <- list_of_mtuples[!vapply(list_of_mtuples, is.null,
                                                 logical(1L))]
  # TODO: do.call(c, list_of_mtuples) ought to work
  Reduce(c, list_of_mtuples)
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
            callNextMethod()
            margin <- "  "
            cat(margin, "methinfo: ", summary(methinfo(object)), "\n", sep = "")
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
            if (missing(x)) {
              args <- unname(list(...))
            } else {
              args <- unname(list(x, ...))
            }
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