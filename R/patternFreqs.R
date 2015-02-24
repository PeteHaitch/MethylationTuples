### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### patternFreqs: Estimate the relative frequencies of methylation patterns from 
### a MethPat object.
###

# TODO: Document
# 
#' A naive estimate of the frequency of methylation patterns.
#' @aliases patternFreqs
#'
#' @export
patternFreqs <- function(methpat, min_cov = 5L, ...) {
  
  # TODO: The requirement of 1 sample allows me to greatly reduce the size of 
  # the intermediate methpat object on which the computations are performed by 
  # removing all m-tuples with coverage < min_cov.
  if (ncol(methpat) > 1) {
    stop(paste0("Only 'MethPat' objects containing data on a single sample ", 
                "are currently supported."))
  }
  # TODO: If stranded MethPat objects are allowed then need to include strand in 
  # returned value.
  if (.stranded(methpat)) {
    stop(paste0("Only 'MethPat' objects that have been processed with ", 
                "'MethylationTuples::collapseStrand' are currently supported"))
  }
  if (size(methpat) < 3) {
    warning(paste0("It is recommended that only 'MethPat' objects with 'size' ", 
                   "> 2 are used for estimating the relative frequencies of ", 
                   "methylation patterns."))
  }
  # TODO: Remove as.vector() if allowing multiple samples
  cov <- as.vector(getCoverage(methpat))
  methpat <- methpat[!is.na(cov) & cov >= min_cov, ]
  cov <- getCoverage(methpat)
  freq <- lapply(assays(methpat), function(assay, cov) {
    assay / cov
  }, cov = cov)
  freq <- matrix(unlist(freq, use.names = FALSE), nrow = length(freq), 
                 byrow = TRUE)
  # Must transpose the resulting column-sorted matrix
  freq <- t(.Call(Cpp_MethylationTuples_sortMatrix, 
                  PACKAGE = "MethylationTuples", 
                  x = freq, sort_direction = "descend", dim = 0L))
  chr <- data.table(chr = as.vector(seqnames(methpat)))
  ipd <- setnames(as.data.table(IPD(methpat)), 
                  paste0("IPD", seq_len(size(methpat) - 1L)))
  freq <- as.data.table(freq)
  setnames(freq, paste0("h", seq_len(ncol(freq))))
  sample <- data.table(sample = colnames(methpat))
  cbind(chr, ipd, sample, freq)
}
