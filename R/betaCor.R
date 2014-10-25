### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### bestCor: Compute within-sample correlations of beta-values.
###

# TODO: Get confidence interval for correlations, e.g., via cor.test
# TODO: Update docs

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
#' @param min_cov An \code{integer} specifying the minimum sequencing coverage 
#' required to use a beta-value.
#' @param pair_type A character string giving the type of pairs to be 
#' constructed when computing correlations. One of "\code{neighbours}" or "\code{all}", can be abbreviated. Please see the below section, "Construction 
#' of pairs of beta-values", for details.
#' @param ref_loci An \code{\link{MTuples}} object with the locations of all 
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
betaCor <- function(methpat, pair_type = c('adjacent', 'all', 'ref_adjacent'), 
                    ipd = seq_len(2000L), ref_loci,
                    method = c('pearson', 'spearman'), 
                    min_cov = 5L,
                    feature) {
  if (!is(methpat, "MethPat") || size(methpat) != 1L) {
    stop("'methpat' must be a 'MethPat' object containing 1-tuples.")
  }
  pair_type <- match.arg(pair_type)
  if (pair_type == 'all') {
    if (!is.numeric(ipd) || !is.vector(ipd)) {
      stop("'ipd' must be an integer vector.")
    }
    ipd <- as.integer(ipd)
  } else if (pair_type == 'adjacent') {
    # NOTHING TODO
  } else if (pair_type == "ref_adjacent") {
    stop("Sorry, 'pair_type = \"ref_adjacent\"' not yet implemented.")
    # TODO: Check how many samples in methpat. It's easiest to implement 
    # ref_adjacent if there is only 1 sample in methpat. Otherwise have to have 
    # a different ref_loci for each sample and deal with that.
    if (missing(ref_loci) || !is(ref_loci, "MTuples") || size(ref_loci) != 1) {
      stop("If 'pair_type' = 'ref_adjacent', then must supply 'ref_loci'.")
      seqinfo <- try(merge(seqinfo(methpat), seqinfo(ref_loci)), silent = TRUE)
      if (is(seqinfo, "try-error")) {
        # TODO: Stricter check of seqinfo compatability, e.g., identical?
        stop("'methpat' and 'ref_loci' have incompatible 'seqinfo'.")
      }
      if (!all(methtype(methpat) %in% methtype(ref_loci))) {
        stop("'methpat' and 'mtuples' have incompatible 'seqinfo'.")
      }
      # TODO: This is very slow. Should it return an error or report the loci 
      # that have no match (or their count) 
      if (isTRUE(any(!overlapsAny(methpat, ref_loci)))) {
        stop("All loci in 'methpat' must also be present in 'ref_loci'.")
      }
    }
  } 
  method <- match.arg(method)
  if (!missing(feature)) {
    if (!is(feature, "GRanges") || !isDisjoint(feature)) {
      stop("'feature' must be a 'GRanges' object with disjoint ranges.")
    }
  }
  
  # Order methpat and extract sorted rowData
  if (pair_type == 'ref_adjacent') {
    # TODO: Add "missing" loci to methpat, i.e., those loci that have 
    # insufficient sequencing coverage in the sample. 
    # E.g., methpat <- rbind(methpat, MethPat(loci_not_in_methpat))
    # TODO: Find out if is.unsorted works on GTuples Could save an expensive 
    # sort if data are already sorted.
    #     if (is.unsorted(ref_loci)) {
    #       ref_loci <- sort(ref_loci)
    #     }
    ref_loci <- sort(ref_loci)
  }
  # TODO: Check whether is.unsorted works on GTuples. Could save an expensive 
  # sort if data are already sorted.
  #   if (is.unsorted(rowData(methpat))) {
  #     methpat_order <- order(methpat)
  #     methpat_rd_sorted <- rowData(methpat)[methpat_order]
  #   } else {
  #     methpat_order <- seq_len(nrow(methpat))
  #     methpat_rd_sorted <- rowData(methpat)
  #   }
  methpat_order <- order(methpat)
  methpat_rd_sorted <- rowData(methpat)[methpat_order]
  betas <- betaVal(methpat, min_cov)
  
  # Get the "feature status" (feature_status) of each loci.
  # feature_status = TRUE if overlaps feature, FALSE otherwise
  if (!missing(feature)) {
    feature_status <- overlapsAny(methpat_rd_sorted, feature)
  } else {
    feature_status <- rep(FALSE, length(methpat_rd_sorted))
  }
  # Create map between IPD-strand-feature_status and an integer ID.
  # Need to define possible IPDs in order to create map.
  if (pair_type == "adjacent" || pair_type == "ref_adjacent") {
    ipd <- sort(unique(diff(start(methpat_rd_sorted))))
    ipd <- ipd[ipd > 0]
  }
  id_dt <- setDT(expand.grid(IPD = ipd, 
                             strand = levels(strand(methpat)),
                             pair_feature_status = 0:3))
  id_dt[, c("KEY", "ID") := list(paste(IPD, strand, pair_feature_status, 
                                       sep = ''),
                                 seq_len(nrow(id_dt)))]
  setkey(id_dt, ID)
  
  # Create pairs of beta-values
  if (pair_type == 'all') {
    # TODO: Benchmark and profile. 
    # 2-3 hours for a MethPat object with 1 sample + 56M CpGs which makes 
    # 1.6 billion pairs and is ~40GB in size.
    pairs <- setkey(setDT(.Call(Cpp_MethylationTuples_makeAllPairs, 
                                methpat_order,
                                as.character(seqnames(methpat_rd_sorted)), 
                                as.character(strand(methpat_rd_sorted)),
                                start(methpat_rd_sorted),                                 
                                feature_status,
                                ipd,
                                betas, 
                                id_dt)), ID, sample)    
  } else if (pair_type == 'adjacent' || pair_type == 'ref_adjacent') {
    # TODO: Benchmark and profile. 
    # 4-5 minutes for a MethPat object with 3 samples and 54 million CpGs.
    # Returns a data.table with 162 million rows and is ~4GB in size.
    pairs <- setkey(setDT(.Call(Cpp_MethylationTuples_makeAdjacentPairs, 
                                methpat_order,
                                as.character(seqnames(methpat_rd_sorted)), 
                                as.character(strand(methpat_rd_sorted)),
                                start(methpat_rd_sorted),                                
                                feature_status,
                                betas,
                                id_dt)), ID, sample)
  }
  
  # Compute correlations
  cors <- pairs[, list(cor = suppressWarnings(cor(beta1, beta2, 
                                                  use = "na.or.complete", 
                                                  method = method))), 
                by = list(ID, sample)]
  
  # Join cors and id_dt. Add sample names back.
  val <- id_dt[cors]
  val[, c("ID", "KEY") := list(NULL, NULL)][, sample := colnames(methpat)[val$sample]]
  if (missing(feature)) {
    val[, pair_feature_status := NULL]
  }
  return(val)
}