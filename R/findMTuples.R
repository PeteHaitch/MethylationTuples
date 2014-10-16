# TODO: Function to find MTuples in a BSgenome

#' @param bsgenome A \code{\link[BSgenome]{BSgenome}} object for the reference 
#' genome of interest.
#' @param methinfo A \code{\link{MethInfo}} object containing the type of 
#' methylation loci for which to search.
#' @param size An \code{integer} specifying the size of tuples for which to 
#' search.
#' @param ... Additional arguments. For example, 
#' you may wish to only search a specific chromosome in the reference genome 
#' via the \code{exclude} argument, which is passed to the 
#' \code{\link[BSgenome]{BSParams}} initializer.
findMTuples <- function(bsgenome, methinfo, size, ...) {
  # 1. Argument checks
  
  # 2. Find all methylation loci with type given by methinfo
  np <- length(methtype(methinfo))
  pattern <- DNAStringSet(methtype(methinfo))
  pattern <- c(pattern, reverseComplement(pattern))
  strand <- rep(c('+', '-'), each = np)
  names(pattern) <- strand
  names(strand) <- strand
  # TODO: Why doesn't bpmapply work?
  # TODO: Use explicit loops and reverse order of loops, i.e. outer loop 
  # chromosomes and inner loop patterns (see http://www.bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/GenomeSearching.pdf)
  loci <- mapply(FUN = function(bsgenome, pattern_, strand_) {
    
    bsparams <- new("BSParams", X = bsgenome, 
                    FUN = function(pattern, subject, strand, ...) {
                      if (strand == '+') {
                        start(matchPattern(pattern, subject, fixed = FALSE, 
                                           ...))
                      } else if (strand == '-') {
                        end(matchPattern(pattern, subject, fixed = FALSE, ...))
                      }
                    }, simplify = TRUE, ...)
    bsapply(bsparams, pattern = pattern_, strand = strand_)
  }, pattern_ = pattern, strand_ = strand, 
  MoreArgs = list(bsgenome = bsgenome), SIMPLIFY = FALSE, ...)
  # Now, combine methylation loci on a per-chromosome + per-strand basis.
  # TODO: First, benchmark to see if it's worth fixing. It seems that there 
  # must be a better way to do this, e.g. pre-allocating vectors or lapply 
  # method. 
  pps <- length(pattern) / 2
  fwd_loci <- loci[[1]]
  rev_loci <- loci[[1 + pps]]
  if (pps > 1) {
    for (i in seq.int(from = 2, to = pps, by = 1)) {
      for (j in seq_along(loci[[i]])) {
        fwd_loci[[j]] <- c(fwd_loci[[j]], loci[[i]][[j]])
      }
    }
    for(i in seq.int(from = pps + 2, to = 2 * pps, by = 1)) {
      for (j in seq_along(loci[[i]])) {
        rev_loci[[j]] <- c(rev_loci[[j]], loci[[i]][[j]])
      }
    }
  }
  fwd_loci <- lapply(fwd_loci, sort)
  rev_loci <- lapply(rev_loci, sort)
  
  # 3. Create all m-tuples of those loci
  # TODO: Is there a better way to create a single MTuples object than 
  # unlist(MTuplesList(list_of_MTuples))?
  mt <- bpmapply(function(seqname, fwd_loci_, rev_loci_, m, methinfo, seqinfo) {
    strand <- c(Rle('+', length(fwd_loci_) - m + 1), 
                Rle('-', length(rev_loci_) - m + 1))
    # TODO: Generators would be cool here
    fwd_idx <- sapply(seq_len(m), function(i, m, n) {
      seq.int(i, n + i - m, 1)
    }, m = m, n = length(fwd_loci_))
    rev_idx <- sapply(seq_len(m), function(i, m, n) {
      seq.int(i, n + i - m, 1)
    }, m = m, n = length(rev_loci_))
    fwd_tuples <- matrix(fwd_loci_[fwd_idx], ncol = m)
    rev_tuples <- matrix(rev_loci_[rev_idx], ncol = m)
    MTuples(GTuples(seqnames = seqname, tuples = rbind(fwd_tuples, rev_tuples), 
                    strand = strand, seqinfo = seqinfo), methinfo = methinfo)
  }, seqname = names(fwd_loci), fwd_loci_ = fwd_loci, 
  rev_loci_ = rev_loci, 
  MoreArgs = list(m = m, methinfo = methinfo, seqinfo = seqinfo(bsgenome)), 
  SIMPLIFY = FALSE)
  unlist(MTuplesList(mt), use.names = FALSE)
}