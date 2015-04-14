### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### methLevelCor: Compute within-sample correlations of methylation levels.
###

# TODO: Get confidence interval for Spearman and Kendall correlation 
# coefficients.
# TODO: Specify or recommend a value of min_cov?
# TODO: Should methLevelCor return unstratified estimates along with stratified 
# estimates? Not hard, just join id_dt[pairs] and then compute correlations 
# ignoring pair_feature_status.
# TODO: Update docs

#' Compute within-sample correlations of pairs of methylation levels.
#' 
#' Given a \code{\link{MethPat}} object containing 1-tuples, 
#' \code{methLevelCor} computes the within-sample correlations of pairs of 
#' methylation levels. Methylation levels are computed by 
#' \code{\link{methLevel}}. 
#' 
#' The manner in which the pairs are constructed is determined by additional 
#' arguments - see the section, "Constructing pairs of methylation loci". 
#' Correlations are stratified by the \code{strand}, the intra-pair distance 
#' and the \code{feature} (if supplied) of the pairs.
#' 
#' @param methpat A \code{\link{MethPat}} object containing 1-tuples.
#' @param pair_type A character string giving the type of pairs to be 
#' constructed when computing correlations. One of "\code{neighbours}" or 
#' "\code{all}", can be abbreviated. Please see the below section, "Construction 
#' of pairs of methylation loci", for details.
#' @param ref_loci An \code{\link{MTuples}} object with the locations of all 
#' methylation loci 1-tuples in the sample/reference genome. Only required if 
#' \code{pair_type = "neighbours"} and otherwise ignored. Please see the below 
#' section, "Construction of pairs of methylation loci", for details.
#' @param max_ipd An \code{integer} specifying the maximal IPD of pairs 
#' (default: 2000). Only required if \code{pair_type = "all"} and otherwise 
#' ignored. Please see the below section, "Construction of pairs of 
#' methylation loci", for details.
#' @param method A character string indicating which correlation coefficient is 
#' to be computed. One of "\code{pearson}" (default) or "\code{spearman}", can 
#' be abbreviated.
#' @param feature An optional \code{\link[GenomicRanges]{GRanges}} object with 
#' the locations of a "genomic feature". This 
#' \code{\link[GenomicRanges]{GRanges}} object must be disjoint (see 
#' \code{\link[GenomicRanges]{isDisjoint}}). Please see the 
#' below section, "Stratifying pairs by a genomic feature", for details.
#' @param ... Additional arguments passed to \code{\link{methLevel}}, such as 
#' \code{min_cov}, \code{statistic} and \code{offset}.
#' 
#' @section Constructing pairs of methylation loci:
#' There are two algorithms for constructing pairs of methylation loci: 
#' \code{neighbours}, which requires the specification of \code{mtuples}, and 
#' \code{all}, which requires the specification of \code{max_ipd}.
#' 
#' \itemize{
#'  \item{\code{pair_type = "neighbours"}:}{ Creates pairs of methylation 
#'  levels from neighbouring methylation loci. It checks that the constructed 
#'  pairs are neighbours by comparing these to the set of all known methylation 
#'  loci in the sample, hence this must be given by \code{mtuples}.}
#'  
#'  \item{\code{pair_type = "all"}:}{ Creates pairs of methylation levels using 
#'  all methylation loci in the sample such that each pair has an intra-pair 
#'  distance less than or equal to \code{max_ipd}.}
#' 
#' }
#' 
#' @section Stratifying pairs by a genomic feature:
#' Pairs of methylation levels may be stratified by whether they are inside 
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
#' make \code{methLevelCor} fairly general, if you have more complex 
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
methLevelCor <- function(methpat, 
                         pair_type = c('adjacent', 'all', 'strict_adjacent'), 
                         ipd = seq_len(2000L), ref_loci,
                         method = c("pearson", "kendall", "spearman"), 
                         conf.level = 0.95,
                         feature, 
                         ...) {
  
  if (nrow(methpat) > .Machine$integer.max) {
    stop(paste0("Sorry, 'methLevelCor' doesn't yet support 'methpat' objects ", 
                "with more than ", .Machine$integer.max, 
                " (.Machine$integer.max) rows."))
  }
  
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
  } else if (pair_type == "strict_adjacent") {
    # TODO: Implement using sample-specific loci.
    # TODO: Check how many samples in methpat. It's easiest to implement 
    # strict_adjacent if there is only 1 sample in methpat. Otherwise have to have 
    # a different ref_loci for each sample and deal with that.
    warning(paste0("'pair_type = \"strict_adjacent\"' is currently based on ",
                   "reference loci, not sample-specific loci."))
    if (missing(ref_loci) || !is(ref_loci, "MTuples") || size(ref_loci) != 1L) {
      stop("If 'pair_type' = 'strict_adjacent', then must supply 'ref_loci'.")
    }
    seqinfo <- try(merge(seqinfo(methpat), seqinfo(ref_loci)), silent = TRUE)
    if (is(seqinfo, "try-error")) {
      # TODO: Stricter check of seqinfo compatability, e.g., identical?
      stop("'methpat' and 'ref_loci' have incompatible 'seqinfo'.")
    }
    if (!all(methtype(methpat) %in% methtype(ref_loci))) {
      stop("'methpat' and 'mtuples' have incompatible 'methinfo'.")
    }
    # TODO: This is very slow. Should it return an error or report the loci 
    # that have no match (or their count) 
    if (isTRUE(any(!overlapsAny(methpat, ref_loci)))) {
      warning("Some loci in 'methpat' not present in 'ref_loci'.")
    }
  }
  method <- match.arg(method)
  if (!missing(feature)) {
    if (!is(feature, "GRanges") || !isDisjoint(feature)) {
      stop("'feature' must be a 'GRanges' object with disjoint ranges.")
    }
  }
  
  # Order methpat and extract sorted rowRanges
  if (pair_type == 'strict_adjacent') {
    # Add "missing" loci to methpat, i.e., those loci that have 
    # insufficient sequencing coverage in the sample.
    missing_loci <- ref_loci[!overlapsAny(ref_loci, methpat, type = "equal")]
    missing_loci <- suppressWarnings(
      MethPat(assays = SimpleList(M = matrix(NA_integer_, 
                                             nrow = length(missing_loci),
                                             ncol = ncol(methpat), 
                                             dimnames = 
                                               list(NULL, colnames(methpat))),
                                  U = matrix(NA_integer_, 
                                             nrow = length(missing_loci),
                                             ncol = ncol(methpat), 
                                             dimnames = 
                                               list(NULL, colnames(methpat)))),
              rowRanges = missing_loci, 
              colData = colData(methpat))
    )
    methpat <- rbind(methpat, missing_loci)
  }
  # TODO: Check whether is.unsorted works on GTuples. Could save an expensive 
  # sort if data are already sorted.
  #   if (is.unsorted(rowRanges(methpat))) {
  #     methpat_order <- order(methpat)
  #     methpat_rd_sorted <- rowRanges(methpat)[methpat_order]
  #   } else {
  #     methpat_order <- seq_len(nrow(methpat))
  #     methpat_rd_sorted <- rowRanges(methpat)
  #   }
  methpat_order <- order(methpat)
  methpat_rd_sorted <- rowRanges(methpat)[methpat_order]
  meth_level <- methLevel(methpat, ...)
  
  if (!missing(feature)) {
    in_feature <- overlapsAny(methpat_rd_sorted, feature)
  } else {
    in_feature <- rep(NA, length(methpat_rd_sorted))
  }
  # Create map between IPD-strand-in_feature and an integer ID.
  # Need to define possible IPDs in order to create map.
  if (pair_type == "adjacent" || pair_type == "strict_adjacent") {
    ipd <- sort(unique(diff(start(methpat_rd_sorted))))
    ipd <- ipd[ipd > 0]
  }
  in_feature_levels <- unique(in_feature)
  pair_feature_status <- sort(unique(rowSums(expand.grid(in_feature_levels, 
                                                         in_feature_levels))),
                              na.last = FALSE)
  id_dt <- setDT(expand.grid(IPD = ipd, 
                             strand = levels(droplevels(strand(methpat))),
                             pair_feature_status = pair_feature_status))
  id_dt[, c("KEY", "ID") := list(paste(IPD, strand, pair_feature_status, 
                                       sep = ''),
                                 seq_len(nrow(id_dt)))]
  setkey(id_dt, ID)
  
  # Create pairs of methylation levels
  if (pair_type == 'all') {
    # TODO: Benchmark and profile.
    # ~2.7 hours for a MethPat object with 54,256,149 CpGs and 3 samples 
    # stratified by CGI-status.
    pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                       methpat_order,
                       as.character(seqnames(methpat_rd_sorted)), 
                       as.character(strand(methpat_rd_sorted)),
                       start(methpat_rd_sorted),                                 
                       in_feature,
                       ipd,
                       id_dt)
  } else if (pair_type == 'adjacent' || pair_type == 'strict_adjacent') {
    # TODO: Benchmark and profile.
    # ~7 minutes for a MethPat object with 54,256,149 CpGs and 3 samples 
    # stratified by CGI-status.
    pairs_idx <- .Call(Cpp_MethylationTuples_makeAdjacentPairs, 
                   methpat_order,
                   as.character(seqnames(methpat_rd_sorted)), 
                   as.character(strand(methpat_rd_sorted)),
                   start(methpat_rd_sorted),                                
                   in_feature,
                   id_dt)
  }
  
  # Compute correlations
  # TODO: bplapply.
  cors_list <- lapply(colnames(methpat), function(sample_name, 
                                                  pairs_idx, meth_level) {
    meth_level_pairs <- data.table(
      ID = pairs_idx[["ID"]], 
      sample = sample_name,
      meth_level_1 = meth_level[pairs_idx[["i"]], sample_name], 
      meth_level_2 = meth_level[pairs_idx[["j"]], sample_name])
    meth_level_pairs[, .myCor(meth_level_1, meth_level_2, method = method, 
                         conf.level = conf.level), by = list(ID, sample)]
  }, pairs_idx = pairs_idx, meth_level = meth_level)
  cors <- setkey(rbindlist(cors_list), ID, sample)
  
  # Join cors and id_dt. Add sample names back.
  # TODO: Check that this is doing the correct sort of join, e.g., inner-join, 
  # outer-join, etc.
  val <- id_dt[cors]
  val[, c("ID", "KEY") := list(NULL, NULL)]
  return(val)
}
