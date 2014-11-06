### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeAdjacentPairsCpp
###

context("makeAdjacentPairsCpp")

test_that("Works with no feature", {
  ipd <- sort(unique(diff(start(methpat_rd_sorted))))
  ipd <- ipd[ipd > 0]
  in_feature <- rep(NA, 13)
  in_feature_levels <- unique(in_feature)
  pair_feature_status <- sort(unique(rowSums(expand.grid(in_feature_levels, 
                                                         in_feature_levels))),
                              na.last = FALSE)
  id_dt <- setDT(expand.grid(IPD = ipd, 
                             strand = levels(strand(methpat)),
                             pair_feature_status = pair_feature_status))
  id_dt[, c("KEY", "ID") := list(paste(IPD, strand, pair_feature_status, 
                                       sep = ''),
                                 seq_len(nrow(id_dt)))]
  setkey(id_dt, ID)
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAdjacentPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                
                     in_feature,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], 
                   c("19+NA", "10+NA", "10-NA", "10-NA", "5+NA", "5+NA", 
                     "10-NA", "10-NA"))
  expect_identical(pairs_idx$i, methpat_order[c(1, 2, 4, 5, 8, 9, 11, 12)])
  expect_identical(pairs_idx$j, methpat_order[c(2, 3, 5, 6, 9, 10, 12, 13)])
  
})

test_that("Works with feature", {
  ipd <- sort(unique(diff(start(methpat_rd_sorted))))
  ipd <- ipd[ipd > 0]
  in_feature <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, 
                  TRUE, FALSE, FALSE, FALSE)
  in_feature_levels <- unique(in_feature)
  pair_feature_status <- sort(unique(rowSums(expand.grid(in_feature_levels, 
                                                         in_feature_levels))),
                              na.last = FALSE)
  id_dt <- setDT(expand.grid(IPD = ipd, 
                             strand = levels(strand(methpat)),
                             pair_feature_status = pair_feature_status))
  id_dt[, c("KEY", "ID") := list(paste(IPD, strand, pair_feature_status, 
                                       sep = ''),
                                 seq_len(nrow(id_dt)))]
  setkey(id_dt, ID)
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAdjacentPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                
                     in_feature,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], 
                   c("19+2", "10+2", "10-1", "10-2", "5+0", "5+1", "10-0", 
                     "10-0"))
  expect_identical(pairs_idx$i, methpat_order[c(1, 2, 4, 5, 8, 9, 11, 12)])
  expect_identical(pairs_idx$j, methpat_order[c(2, 3, 5, 6, 9, 10, 12, 13)])
})

# TODO: Tests of methLevelCor functionality aside from "makePairs" functions.