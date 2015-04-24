### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mantelhaen: Estimate within-sample, within-fragment co-methylation.
###

# WARNING: This is *really* bare bones and should be used with caution. 
# It will eventually be polished and made part of the cometh() function.
# TODO: Merge cometh() and mantelhaen(), i.e., make mantelhaen() an option 
# in cometh() via the 'method' argument.
#' @export
mantelhaen <- function(methpat, key = c("seqnames", "IPD"), K = NULL, 
                       BPPARAM = bpparam()) {
  
  stopifnot(size(methpat) == 2L)
  
  mh <- function(MM, MU, UM, UU) {
    x <- array(rbind(MM, MU, UM, UU), dim = c(2, 2, length(MM)))
    val <- mantelhaen.test(x)
    list(estimate = val$estimate,
         CI_lower = val$conf.int[1],
         CI_upper = val$conf.int[2]) 
  }
  
  bplapply(seq_len(ncol(methpat)), function(i, methpat, key, K) {
    assays <- assays(methpat, withDimnames = FALSE)
    y <- data.table(seqnames = as.factor(seqnames(methpat)),
                    IPD = as.vector(IPD(methpat)),
                    MM = assays[["MM"]][, i, drop = TRUE],
                    MU = assays[["MU"]][, i, drop = TRUE],
                    UM = assays[["UM"]][, i, drop = TRUE],
                    UU = assays[["UU"]][, i, drop = TRUE],
                    key = key)
    y <- y[!is.na(MM) & (MM + MU + UM + UU) > 1, ]
    y <- y[, n := .N, by = key(y)][n > 1, ]
    if (!is.null(K)) {
      y[, bin := rep(seq_len(n), each = K), by = key]
      setkeyv(y, c(key, "bin"))
      y <- y[, n := .N, by = key(y)][n > 1, ]
    }
    
    y[, mh(MM, MU, UM, UU), by = key(y)]
  }, methpat = methpat, key = key, K = K, BPPARAM = BPPARAM)
}