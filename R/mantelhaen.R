### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mantelhaen: Estimate within-sample, within-fragment co-methylation.
###

# Mantel-Haenszel estimator of a common log odds ratio of a 2x2xK table
# Also estimates 95% CI.
MH <- function(x, offset = 0) {
  x <- x + offset
  MM <- x[2, 2, ]
  MU <- x[2, 1, ]
  UM <- x[1, 2, ]
  UU <- x[1, 1, ]
  d <- apply(x, 3, sum)
  or <- sum((UU * MM) / d) / sum((MU * UM) / d)
  # One or more margins must be identically zero so add offset
  if (is.nan(or) || isTRUE(all.equal(or, 0)) || is.infinite(or)) {
    warning(paste0("One or more cell is (probably) identically zero in all K ", 
                   "2x2 tables. Adding offset = 0.5"))
    MH(x, 0.5)
  } else {
    log_ase2 <- sum((UU + MM) * (UU * MM) / (d ^ 2)) / 
      (2 * (sum(UU * MM / d)) ^ 2) +
      sum(((UU + MM) * (MU * UM) + (UM + MU) * (UU * MM)) / d ^ 2) /
      (2 * sum(UU * MM / d) * sum(MU * UM / d)) +
      sum((UM + MU) * (UM * MU) / (d ^ 2)) / (2 * (sum(UM * MU / d)) ^ 2)
    ci <- c(exp(log(or) - qnorm(0.975) * sqrt(log_ase2)), 
            exp(log(or) + qnorm(0.975) * sqrt(log_ase2)))
    list(estimate = or, conf.int = ci)
  }
}


# WARNING: This is *really* bare bones and should be used with caution. 
# It will eventually be polished and made part of the cometh() function.
# TODO: Merge cometh() and mantelhaen(), i.e., make mantelhaen() an option 
# in cometh() via the 'method' argument.
#' @export
mantelhaen <- function(methpat, key = c("seqnames", "IPD"), K = NULL, 
                       BPPARAM = bpparam()) {
  
  stopifnot(size(methpat) == 2L)

  
  bplapply(seq_len(ncol(methpat)), function(i, methpat, key, K) {
    
    mh <- function(MM, MU, UM, UU) {
      x <- array(rbind(MM, MU, UM, UU), dim = c(2, 2, length(MM)))
      val <- MH(x)
      list(estimate = val$estimate,
           CI_lower = val$conf.int[1],
           CI_upper = val$conf.int[2]) 
    }
    
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