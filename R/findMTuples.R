# TODO: Should the example use library or require?
# TODO: Unit tests. This could be difficult since the function requires a 
# large input file (bsgenome) and takes a while to run; what's best practice?

#' Find tuples of methylation loci of a given size in a reference genome.
#' 
#' @param bsgenome A \code{\link[BSgenome]{BSgenome}} object for the reference 
#' genome of interest.
#' @param methinfo A \code{\link{MethInfo}} object containing the type of 
#' methylation loci for which to search.
#' @param size An \code{integer} specifying the size of tuples for which to 
#' search.
#' @param exclude A character vector of \code{seqnames} to be filtered out from 
#' the \code{bsgenome}.
#' 
#' @return An \code{\link{MTuples}} object containing the tuples.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Find all 2-tuples of CG in the C. elegans reference genome
#' library(BSgenome.Celegans.UCSC.ce2)
#' findMTuples(bsgenome = Celegans, methinfo = MethInfo('CG'), size = 2)
#' }
findMTuples <- function(bsgenome, methinfo, size, exclude = NULL) {
  # 0. Argument checks
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
  if (!is.character(exclude) && !is.null(exclude)) {
    stop("'exclude' must be a character vector.")
  }
  
  # 1. Variable initialisation
  methtype <- methtype(methinfo)
  seqinfo <- seqinfo(bsgenome)
  seqnames <- seqnames(bsgenome)
  seqnames <- seqnames[!seqnames %in% exclude]
  list_of_mtuples <- vector(mode = "list", length = length(seqnames))
  names(list_of_mtuples) <- seqnames
  
  # 2. Find all methylation loci with type given by methinfo
  # TODO: foreach would give instant parallelisation; can this be done via
  # BiocParallel and is it worth it?
  for (seqname in seqnames) {
    subject <- bsgenome[[seqname]]
    fwd_loci <- c()
    rev_loci <- c()
    for (mt in methtype) {
      mt <- DNAString(mt)
      rcmt <- reverseComplement(mt)
      tmp_fwd_loci <- start(matchPattern(mt, subject, fixed = FALSE))
      tmp_rev_loci <- end(matchPattern(rcmt, subject, fixed = FALSE))
      fwd_loci <- c(fwd_loci, tmp_fwd_loci)
      rev_loci <- c(rev_loci, tmp_rev_loci)
    }
    strand <- c(Rle('+', length(fwd_loci) - size + 1), 
                Rle('-', length(rev_loci) - size + 1))
    fwd_loci <- sort(fwd_loci)
    rev_loci <- sort(rev_loci)
    # TODO: Iterators/generators would be cool here rather than actually 
    # creating the index variables.
    fwd_idx <- sapply(seq_len(size), function(i, size, n) {
      seq.int(i, n + i - size, 1)
    }, size = size, n = length(fwd_loci))
    rev_idx <- sapply(seq_len(size), function(i, size, n) {
      seq.int(i, n + i - size, 1)
    }, size = size, n = length(rev_loci))
    fwd_tuples <- matrix(fwd_loci[fwd_idx], ncol = size)
    rev_tuples <- matrix(rev_loci[rev_idx], ncol = size)
    tuples <- rbind(fwd_tuples, rev_tuples)
    list_of_mtuples[seqname] <- MTuples(gtuples = GTuples(seqnames = seqname, 
                                                          tuples = tuples, 
                                                          strand = strand, 
                                                          seqinfo = seqinfo),
                                        methinfo = methinfo)
  }
  # TODO: Is there a better way to create a single MTuples object than 
  # unlist(MTuplesList(list_of_MTuples))?
  unlist(MTuplesList(list_of_mtuples), use.names = FALSE)
}