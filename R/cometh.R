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
#' @param ref_loci An \code{\link{MTuples}} object with the locations of all 
#' methylation loci 1-tuples in the sample/reference genome. Only required if 
#' \code{pair_type = "strict_adjacent"} and otherwise ignored. Please see the below 
#' section, "Choice of 'pair_type'", for details.
#' @param method A character string giving the method used to estimate 
#' co-methylation. One of "\code{lor}" or "\code{pearson}", can be abbreviated. 
#' Please see the below section, "Choice of statistic", for details.
#' @param min_cov
#' 
#' 
#' @section Choice of 'pair_type':
#' 
#' @section Choice of statistic:
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
                   offset = 1L) {
  
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
    pair_feature_status <- overlapsAny(methpat, feature, type = 'start') + 
      overlapsAny(methpat, feature, type = 'end')
  } else {
    pair_feature_status <- rep(NA, nrow(methpat))
  }
  
  cov <- getCoverage(methpat)
  if (method == "lor") {
    statistic <- log2((assay(methpat, "MM") * assay(methpat, "UU") + offset) / 
                        (assay(methpat, "MU") * assay(methpat, "UM") + offset))
    # Compute confidence interval
    sigma <- sqrt(1 / (assay(methpat, "MM") + offset) + 
                    1 / (assay(methpat, "MU") + offset) + 
                    1 / (assay(methpat, "UM") + offset) +
                    1 / (assay(methpat, "UU") + offset))
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
    num <- assay(methpat, "MM") * assay(methpat, "UU") - 
      assay(methpat, "UM") * assay(methpat, "MU")
    denom <- sqrt(assay(methpat, "MM") + assay(methpat, "MU")) *
      sqrt((assay(methpat, "UM") + assay(methpat, "UU"))) *
      sqrt((assay(methpat, "MM") + assay(methpat, "UM"))) * 
      sqrt((assay(methpat, "MU") + assay(methpat, "UU")))
    statistic <- num / denom
    # Compute confidence interval using the Fisher transformation
    # TODO: Can produce "In atanh(statistic) : NaNs produced". 
    # Suspect this is due to NA in statistic but haven't been able to 
    # reproduce.
    z <- atanh(statistic)
    sigma <- suppressWarnings(1 / sqrt(cov - 3L))
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
  statistic[cov < min_cov] <- NA
  CI_lower[cov < min_cov] <- NA
  CI_upper[cov < min_cov] <- NA
  
  data.table(chr = as.character(seqnames(methpat)), 
             pos1 = as.integer(start(methpat)), 
             pos2 = as.integer(end(methpat)), 
             strand = as.character(strand(methpat)), 
             pair_feature_status = pair_feature_status, 
             sample = rep(colnames(methpat), each = nrow(methpat)), 
             statistic = as.vector(statistic), 
             CI_lower = as.vector(CI_lower), 
             CI_upper = as.vector(CI_upper))
}
