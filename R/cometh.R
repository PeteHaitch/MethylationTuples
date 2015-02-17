### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cometh: Compute within-sample, within-fragment co-methylation.
###

# TODO: Should pair_type even by an option? That is, should the user have already 
# decided whether to remove pairs of loci with NIL > 0? If so, provide a function 
# filter_nil().
# TODO: All filter_* functions should identically either retain "filtered" or 
# remove "filtered" objects. Use the same logic as dplyr.
# TODO: Finish documenting.

#' Estimate within-sample, within-fragment co-methylation using 2-tuples.
#' 
#' @param methpat A \code{\link{MethPat}} object containing 2-tuples.
#' @param pair_type A character string giving the type of pairs to be 
#' constructed when computing correlations. One of "\code{neighbours}" or 
#' "\code{all}", can be abbreviated. Please see the below section, "Choice of 
#' 'pair_type'", for details.
#' @param method A character string indicating which statistic to use to 
#' estimate co-methylation. This must be (an abbreviation of) one of the 
#' strings "\code{lor}" (\eqn{log_{2}(odds-ratio)}) or 
#' "\code{pearson}" (Pearson product-moment correlation coefficient). Please 
#' see the below section, "Choice of statistic", for details.
#' @param min_cov An \code{integer} specifying the minimum coverage required 
#' in order to use a 2-tuple when estimating co-methylation.
#' @param alternative Indicates the alternative hypothesis and must be one of 
#' "two.sided", "greater" or "less". You can specify just the initial letter. 
#' "greater" corresponds to positive association, "less" to negative 
#' association.
#' @param feature An optional \code{\link[GenomicRanges]{GRanges}} object with 
#' the locations of a "genomic feature". This 
#' \code{\link[GenomicRanges]{GRanges}} object must be disjoint (see 
#' \code{\link[GenomicRanges]{isDisjoint}}). Please see the 
#' below section, "Stratifying pairs by a genomic feature", for details.
#' @param offset A \code{numeric} vector with length 1 used when computing 
#' M-values (default: 0.5).
#' 
#' @section Choice of 'pair_type':
#' 
#' @section Choice of statistic:
#' 
#' @section Stratifying pairs by a genomic feature:
#' 
#' @export
cometh <- function(methpat, 
                   pair_type = c('all', 'strict_adjacent'), 
                   ref_loci,
                   method = c("lor", "pearson"),
                   min_cov = 5L,
                   alternative = c("two.sided", "less", "greater"),
                   conf.level = 0.95,
                   feature, 
                   offset = 0.5) {
  
  # TODO: Check if this likely to occur in practice
  if (nrow(methpat) > .Machine$integer.max) {
    stop(paste0("Sorry, 'methLevelCor' doesn't yet support 'methpat' objects ", 
                "with more than ", .Machine$integer.max, 
                " (.Machine$integer.max) rows."))
  }
  if (!is(methpat, "MethPat") || size(methpat) != 2L) {
    stop("'methpat' must be a 'MethPat' object containing 2-tuples.")
  }
  
  pair_type <- match.arg(pair_type)
  if (pair_type == "strict_adjacent") {
    stop("Sorry, 'pair_type = \"strict_adjacent\"' not yet implemented.")
    # TODO: Check how many samples in methpat. It's easiest to implement 
    # strict_adjacent if there is only 1 sample in methpat. Otherwise have to have 
    # a different ref_loci for each sample and deal with that.
    if (missing(ref_loci)) {
      stop("If 'pair_type' = 'strict_adjacent', then must supply 'ref_loci',") 
    }
    if(!is(ref_loci, "MTuples") || size(ref_loci) != 1) {
      stop(paste0("'ref_loci' must be an 'MTuples' object containing ", 
                  "1-tuples of methylation loci."))
    }
    seqinfo <- try(merge(seqinfo(methpat), seqinfo(ref_loci)), silent = TRUE)
    if (is(seqinfo, "try-error")) {
      # TODO: Stricter check of seqinfo compatability, e.g., identical?
      stop("'methpat' and 'ref_loci' have incompatible 'seqinfo'.")
    }
    if (!all(methtype(methpat) %in% methtype(ref_loci))) {
      stop("'methpat' and 'mtuples' have incompatible 'seqinfo'.")
    }
    # TODO: Check that all loci in methpat are also in ref_loci.
  }
  
  method <- match.arg(method)
  
  min_cov <- as.integer(min_cov)
  
  alternative <- match.arg(alternative)
  
  if (!missing(conf.level) && 
        (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || 
           conf.level > 1)) {
    stop("'conf.level' must be a single number between 0 and 1")
  }
  
  if (!missing(feature)) {
    rd_start <- GRanges(seqnames(methpat), IRanges(start(methpat), width = 1L), 
                        strand(methpat), seqinfo = seqinfo(methpat))
    rd_end <- GRanges(seqnames(methpat), IRanges(end(methpat), width = 1L), 
                      strand(methpat), seqinfo = seqinfo(methpat))
    pair_feature_status <- Rle(overlapsAny(rd_start, feature) + 
                                 overlapsAny(rd_end, feature))
    rm(rd_start, rd_end)
  } else {
    pair_feature_status <- Rle(NA, nrow(methpat))
  }
  
  cov <- getCoverage(methpat)
  # mc = those loci with the minimum coverage. Only compute statistic for these
  # loci.
  mc <- cov >= min_cov & !is.na(cov)
  # TODO: Does a forced rm() actually help or will it be gc-ed anyway?
  rm(cov)
  if (method == "lor") {
    statistic <- log2(assay(methpat, "MM")[mc] + offset) + 
      log2(assay(methpat, "UU")[mc] + offset) -
      log2(assay(methpat, "MU")[mc] + offset) - 
      log2(assay(methpat, "UM")[mc] + offset)
    # Compute confidence interval
    sigma <- sqrt(1 / (assay(methpat, "MM")[mc] + offset) + 
                    1 / (assay(methpat, "MU")[mc] + offset) + 
                    1 / (assay(methpat, "UM")[mc] + offset) +
                    1 / (assay(methpat, "UU")[mc] + offset))
    CI_lower <- switch(alternative, 
                       less = -Inf, 
                       greater = statistic - sigma * qnorm(conf.level), 
                       two.sided = statistic - sigma * 
                         qnorm((1 + conf.level) / 2))
    CI_upper <- switch(alternative, 
                       less = statistic + sigma * qnorm(conf.level), 
                       greater = Inf,
                       two.sided = statistic + sigma * 
                         qnorm((1 + conf.level) / 2))
  } else if (method == "pearson") {
    # Compute Pearson correlation (equivalent to phi coefficient).
    num <- assay(methpat, "MM")[mc] * assay(methpat, "UU")[mc] - 
      assay(methpat, "UM")[mc] * assay(methpat, "MU")[mc]
    denom <- sqrt(assay(methpat, "MM")[mc] + assay(methpat, "MU")[mc]) *
      sqrt((assay(methpat, "UM")[mc] + assay(methpat, "UU")[mc])) *
      sqrt((assay(methpat, "MM")[mc] + assay(methpat, "UM")[mc])) * 
      sqrt((assay(methpat, "MU")[mc] + assay(methpat, "UU")[mc]))
    statistic <- num / denom
    # Compute confidence interval using the Fisher transformation
    # TODO: Can produce "In atanh(statistic) : NaNs produced". 
    # Suspect this is due to NA in statistic but haven't been able to 
    # reproduce.
    z <- atanh(statistic)
    sigma <- suppressWarnings(1 / sqrt(assay(methpat, "MM")[mc] + 
                                         assay(methpat, "MU")[mc] + 
                                         assay(methpat, "UM")[mc] + 
                                         assay(methpat, "UU")[mc] - 3L))
    sigma[is.nan(sigma)] <- NA
    CI_lower <- switch(alternative, 
                       less = -Inf, 
                       greater = z - sigma * qnorm(conf.level), 
                       two.sided = z - sigma * qnorm((1 + conf.level) / 2))
    CI_lower <- tanh(CI_lower)
    CI_upper <- switch(alternative, 
                       less = z + sigma * qnorm(conf.level), 
                       greater = Inf, 
                       two.sided = z + sigma * qnorm((1 + conf.level) / 2))
    CI_upper <- tanh(CI_upper)
  }
  
  # This clunky construction of a data.table is necessary because working with 
  # Rle objects leads to integer overflow for large methpat objects.
  data.table(
    chr = unlist(lapply(seq_len(ncol(mc)), function(j, mc, sn) {
      as.character(sn[mc[, j]])
    }, mc = mc, sn = seqnames(methpat)), use.names = FALSE), 
    pos1 = unlist(lapply(seq_len(ncol(mc)), function(j, mc, s) {
      s[mc[, j]]
    }, mc = mc, s = start(methpat)), use.names = FALSE), 
    pos2 = unlist(lapply(seq_len(ncol(mc)), function(j, mc, e) {
      e[mc[, j]]
    }, mc = mc, e = end(methpat)), use.names = FALSE), 
    strand = unlist(lapply(seq_len(ncol(mc)), function(j, mc, strand) {
      as.character(strand[mc[, j]])
    }, mc = mc, strand = strand(methpat)), use.names = FALSE),
    pair_feature_status = unlist(lapply(seq_len(ncol(mc)), 
                                        function(j, mc, pfs) {
                                          as.vector(pfs[mc[, j]])
                                        }, mc = mc, pfs = pair_feature_status), 
                                 use.names = FALSE), 
    sample = unlist(mapply(function(j, sn, mc) {
      rep(sn, sum(mc[, j]))
    }, j = seq_len(ncol(mc)), sn = colnames(methpat), 
    MoreArgs = list(mc = mc)), use.names = FALSE), 
    statistic = statistic, 
    sigma = sigma,
    CI_lower = CI_lower, 
    CI_upper = CI_upper)
}
