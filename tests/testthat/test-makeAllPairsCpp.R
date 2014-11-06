### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeAllPairsCpp
###

context("makeAllPairsCpp")

test_that("Works with no feature", {
  ipd <- 1:30
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], 
                   c("19+NA", "29+NA", "10+NA", "10-NA", "20-NA", "10-NA", 
                     "5+NA", "10+NA", "5+NA", "10-NA", "20-NA", "10-NA"))
  expect_identical(pairs_idx$i, 
                   methpat_order[c(1, 1, 2, 4, 4, 5, 8, 8, 9, 11, 11, 12)])
  expect_identical(pairs_idx$j, 
                   methpat_order[c(2, 3, 3, 5, 6, 6, 9, 10, 10, 12, 13, 13)])
  ipd <- 1:5
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], c("5+NA", "5+NA"))
  expect_identical(pairs_idx$i, 
                   methpat_order[c(8, 9)])
  expect_identical(pairs_idx$j, 
                   methpat_order[c(9, 10)])
  ipd <- 6
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], character(0))
  expect_identical(pairs_idx$i, integer(0))
  expect_identical(pairs_idx$j, integer(0))
})

test_that("Works with feature", {
  ipd <- 1:30
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], 
                   c("19+2", "29+2", "10+2", "10-1", "20-1", "10-2", 
                     "5+0", "10+1", "5+1", "10-0", "20-0", "10-0"))
  expect_identical(pairs_idx$i, 
                   methpat_order[c(1, 1, 2, 4, 4, 5, 8, 8, 9, 11, 11, 12)])
  expect_identical(pairs_idx$j, 
                   methpat_order[c(2, 3, 3, 5, 6, 6, 9, 10, 10, 12, 13, 13)])
  ipd <- 1:5
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], c("5+0", "5+1"))
  expect_identical(pairs_idx$i, methpat_order[c(8, 9)])
  expect_identical(pairs_idx$j, methpat_order[c(9, 10)])
  ipd <- 6
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
  pairs_idx <- .Call(Cpp_MethylationTuples_makeAllPairs, 
                     methpat_order,
                     as.character(seqnames(methpat_rd_sorted)), 
                     as.character(strand(methpat_rd_sorted)),
                     start(methpat_rd_sorted),                                 
                     in_feature,
                     ipd,
                     id_dt)
  expect_identical(id_dt[pairs_idx$ID][, KEY], character(0))
  expect_identical(pairs_idx$i, integer(0))
  expect_identical(pairs_idx$j, integer(0))
})
  
  