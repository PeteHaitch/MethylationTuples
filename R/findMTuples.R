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
  pattern <- c(pattern, complement(pattern))
  # TODO: Check bsgenome is constant across iterations of bpmapply
  loci <- bpmapply(FUN = function(bsgenome, pattern_, strand_) {

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
  }, bsgenome = bsgenome, pattern_ = pattern, 
  strand_ = rep(c('+', '_'), each = np), ...)
  # Now, combine methylation loci on a per-chromosome + per-strand basis.
  # TODO: Take care of strand
  # TODO: Sort result
  mloci <- loci[[1]]
  for (i in seq_along(loci)[-1]) {
    for (j in seq_along(loci[[i]])) {
      mloci[[j]] <- c(mloci[[j]], loci[[i]][[j]])
    }
  }
  
  # 3. Create all m-tuples of those loci
  bpmapply(function(mli, seqname, m, methtype) {
    tuples <- matrix()
    # TODO: strand!
    MTuples(seqnames = seqname, 
            tuples = tuples, 
            strand = strand, 
    )
  }, mli = mloci, seqname = names(mloci), m = m, methtype)
  # TODO: Add seqinfo(bsgenome) to final object 
}