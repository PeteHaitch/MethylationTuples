### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bestCor: Compute within-sample correlations of beta-values.
###

# TODO: Note that care should be taken with strand. Perhaps I should force
# unstrand the feature.
# TODO: Currently I drop pairs that span the boundary of a feature. Should I 
# return these as their own category?

#' Compute within-sample correlations of pairs of beta-values.
#' 
#' Given a \code{\link{MethPat}} object containing 1-tuples, \code{betaCor} 
#' computes the within-sample correlations of pairs of beta-values. The manner 
#' in which the pairs are constructed is determined by additional arguments - 
#' see the section, "Construction of pairs of beta-values". Correlations are 
#' stratified by the \code{strand}, the intra-pair distance and the 
#' \code{feature} (if supplied) of the pairs.
#' 
#' @param methpat A \code{\link{MethPat}} object containing 1-tuples.
#' @param pair_type A character string giving the type of pairs to be 
#' constructed when computing correlations. One of "\code{neighbours}" or "\code{all}", can be abbreviated. Please see the below section, "Construction 
#' of pairs of beta-values", for details.
#' @param mtuples An \code{\link{MTuples}} object with the locations of all 
#' methylation loci 1-tuples in the sample/reference genome. Only required if 
#' \code{pair_type = "neighbours"} and otherwise ignored. Please see the below 
#' section, "Construction of pairs of beta-values", for details.
#' @param max_ipd An \code{integer} specifying the maximal IPD of pairs 
#' (default: 2000). Only required if \code{pair_type = "all"} and otherwise 
#' ignored. Please see the below section, "Construction of pairs of 
#' beta-values", for details.
#' @param method A character string indicating which correlation coefficient is 
#' to be computed. One of "\code{pearson}" (default) or "\code{spearman}", can 
#' be abbreviated.
#' @param feature An optional \code{\link[GenomicRanges]{GRanges}} object with 
#' the locations of a "genomic feature". This 
#' \code{\link[GenomicRanges]{GRanges}} object must be disjoint (see 
#' \code{\link[GenomicRanges]{isDisjoint}}). The \code{feature}Please see the 
#' below section, "Stratifying pairs by a genomic feature", for details.
#' @param same_feature A \code{logical(1)}. When stratifying by \code{feature}, 
#' should a pair be required to be in the "same element" of the feature. Please 
#' see the below section, "Stratifying pairs by a genomic feature", for details.
#' @param feature_name A character string with the name of the feature, which 
#' is used in the output, e.g., "CpG island".
#' 
#' @section Constructing pairs of beta-values:
#' There are two algorithms for constructing pairs of beta-values: 
#' \code{neighbours}, which requires the specification of \code{mtuples}, and 
#' \code{all}, which requires the specification of \code{max_ipd}.
#' 
#' \itemize{
#'  \item{\code{pair_type = "neighbours"}:}{ Creates pairs of beta-values 
#'  from neighbouring methylation loci. It checks that the constructed pairs 
#'  are neighbours by comparing these to the set of all known methylation loci 
#'  in the sample, hence this must be given by \code{mtuples}.}
#'  
#'  \item{\code{pair_type = "all"}:}{ Creates pairs of beta-values using all 
#'  methylation loci in the sample such that each pair has an intra-pair 
#'  distance less than or equal to \code{max_ipd}.}
#' 
#' }
#' 
#' @section Stratifying pairs by a genomic feature:
#' Pairs of beta-values may be stratified by whether they are inside 
#' or outside of a "genomic feature". For our example, we will use CpG islands 
#' as our feature. In this case \code{feature} should be a 
#' \code{\link[GenomicRanges]{GRanges}} object containing all CpG islands in 
#' the genome. 
#' 
#' The assignment of pairs as "inside" or "outside" is suprisingly nuanced. The 
#' \code{same_feature} argument controls how a pair is assigned 
#' as being "inside" or "outside" the feature. 
#'  
#' Please think carefully about what defines a pair as being "inside" and 
#' "outside" a feature for your particular needs. While I have attempted to 
#' make \code{betaCor} fairly general, if you have more complex 
#' requirements then you may find the available options to be insufficient.
#' 
#' \itemize{
#'  \item{\code{same_feature = FALSE}: }{A pair is considered to be "inside" 
#'  the feature if and only if both loci that make up in the pair overlap an 
#'  element of \code{feature}. However, the two loci in each pair may overlap 
#'  different elements of \code{feature}, e.g., different CpG islands. 
#'  Similarly, a pair is considered to be "outside" the feature if and only if
#'  both loci that make up the pair overlap a "gap" between elements of 
#'  \code{feature}. However, the two loci in each pair may overlap different 
#'  "gaps", non-CpG island regions of the genome. This is the default.
#'  }
#'  \item{\code{same_feature = TRUE}: A pair is considered to be "inside" the 
#'  feature only if and only if both loci that make up the pair are in the same 
#'  feature, e.g., the same CpG island. Similarly, a pair is considered 
#'  "outside" the feature only if and only if both loci that make up the pair 
#'  are in the same "gap" between elements of \code{feature}. All other pairs 
#'   are discarded, e.g., those where both loci that make up the pair are in a 
#'   CpG island but lie in distinct CpG islands.}
#' }
#' 
#' @export
betaCor <- function(methpat, pair_type = c('neighbours', 'all'), mtuples,
                    max_ipd = 2000, method = c('pearson', 'spearman'),
                    feature, same_feature = FALSE, feature_name) {
  if (!is(methpat, "MethPat") || size(methpat) != 1L) {
    stop("'methpat' must be a 'MethPat' object containing 1-tuples.")
  }
  pair_type <- match.arg(pair_type)
  if (pair_type == "neighbours") {
    if (missing(mtuples) || !is(mtuples, "MTuples") || size(mtuples) != 1) {
      stop("If 'pair_type' = 'neighbours', then must supply 'mtuples'.")
      seqinfo <-try(merge(seqinfo(methpat), seqinfo(mtuples)), silent = TRUE)
      if (is(seqinfo, "try-error")) {
        stop("'methpat' and 'mtuples' have incompatible 'seqinfo'.")
      }
      if (!all(methtype(methpat) %in% methpat(mtuples))) {
        stop("'methpat' and 'mtuples' have incompatible 'seqinfo'.")
      }
      if (isTRUE(any(countOverlaps(methpat, mtuples) == 0L))) {
        stop("All loci in 'methpat' must also be present in 'mtuples'.")
      }
    }
  } else if (pair_type == 'all') {
    if (!is.numeric(max_ipd)) {
      stop("'max_ipd' must be an integer.")
    }
    max_ipd <- as.integer(max_ipd)
  }
  method <- match.arg(method)
  if (!missing(feature)) {
    if (!is(feature, "GRanges") || !isDisjoint(feature)) {
      stop("'feature' must be a 'GRanges' object with disjoint ranges.")
    }
    if (!isTRUEorFALSE(same_feature)) {
      stop("'same_feature' must be TRUE or FALSE.")
    }
    if (missing(feature_name)) {
      warning(paste0("It is recommended that you supply a 'feature_name'.\n",
                     "Instead, using the default, 'feature'."))
    }
    if (!isTRUEorFALSE(ignore_strand)) {
      stop("'ignore_strand' must be TRUE or FALSE")
    }
  }
  
  # Order methpat and extract sorted rowData
  methpat_rd_order <- order(methpat)
  methpat_rd <- rowData(methpat)[methpat_rd_order]
  mtuples <- sort(mtuples)
  beta <- betaVal(methpat)
  
  # Get indices of pairs
  # Should return a list, l(x, y), where x is index of first loci in pair and
  # y is index of second loci in pair wrt to methpath_rd_order
  if (pair_type == 'neighbours') {
    idx <- .makeNeighbourPairsIdx(methpat_rd, mtuples)
  } else if (pair_type == 'all') {
    idx <- .makeAllPairsIdx(methpat_rd, max_ipd)
  }
  # Create pairs
  pairs <- GRanges(seqnames = seqnames(methpat_rd)[idx$x], 
                   ranges = IRanges(start(methpat_rd)[idx$x],
                                    start(methpat_rd)[idx$y]),
                   strand = strand(methpat_rd)[idx$x])
  
  # TODO: Extend feature_idx to include "first in, last out" and 
  # "first out, last in". Should be doable with compare,IRanges-method
  # [Optional] Stratify pairs by feature
  # feature_idx = 2 (inside), 1 (spanning), 0 (outside)
  if (!missing(feature)) {
    if (same_feature) {
      feature_idx <- overlapsAny(pairs, feature, type = 'within') + 
        overlapsAny(pairs, feature, type = 'any')
    } else {
      feature_idx <- overlapsAny(resize(pairs, width = 1, fix = 'start'),
                                 feature, type = 'within') +
        overlapsAny(resize(pairs, width = 1, fix = 'end'), 
                    feature, type = 'within')
    }
  }
  
  # Compute correlations stratified by IPD and feature_idx
  # Create data.table of IPD, feature, beta-values
  beta_pairs <- data.table(ipd = width(pairs) - 1L, feature = feature_idx, 
                           beta_1 = beta[idx$x], beta_2 = beta[idx$y])
  beta_pairs[feature != 1][, list(cor = cor(beta_1, beta_2)), by = list(ipd, feature)]
}