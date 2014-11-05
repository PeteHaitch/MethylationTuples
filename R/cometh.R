### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### cometh: Compute within-sample, within-fragment co-methylation.
###

# TODO: Check method argument
# TODO: Implement method = 'pearson', which computes the correlation Pearson 
# correlation coefficient of each 2x2 table cf. Landan et al. NB: The Pearson 
# correlation coefficient is equivalent to computing the phi coefficient of the 
# 2x2 table.

#' @export
cometh <- function(methpat, 
                   pair_type = c('ref_adjacent', 'all'), 
                   ref_loci,
                   method = c("lor", "pearson"),
                   min_cov = 5L,
                   conf.level = 0.95,
                   feature, 
                   offset = 1) {
  
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
  if (pair_type == "ref_adjacent") {
    stop("Sorry, 'pair_type = \"ref_adjacent\"' not yet implemented.")
    # TODO: Check how many samples in methpat. It's easiest to implement 
    # ref_adjacent if there is only 1 sample in methpat. Otherwise have to have 
    # a different ref_loci for each sample and deal with that.
    if (missing(ref_loci)) {
      stop("If 'pair_type' = 'ref_adjacent', then must supply 'ref_loci',") 
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
  
  if (!missing(feature)) {
    pair_feature_status <- overlapsAny(start(methpat), feature) + 
      overlapsAny(end(methpat), feature)
  } else {
    pair_feature_status <- rep(NA, nrow(methpat))
  }
    
  lor <- log2((assay(methpat, "MM") * assay(methpat, "UU") + offset) / 
                (assay(methpat, "MU") * assay(methpat, "UM") + offset))
  lor[getCoverage(methpat) < min_cov] <- NA
  
  # TODO: Add CI_lower and CI_upper
  data.table(IPD = as.vector(IPD(methpat)), 
             strand = as.character(strand(methpat)), 
             pair_feature_status = pair_feature_status, 
             sample = rep(colnames(methpat), each = nrow(methpat)), 
             lor = as.vector(lor), CI_lower = NA, CI_upper = NA)  
}