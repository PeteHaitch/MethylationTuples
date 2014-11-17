# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("MethPat validity methods")

test_that(".valid.MethPat.rowData works for empty MTuples", {
  expect_error(MethPat(rowData = granges(gt0)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'MTuples' object."))
})

test_that(".valid.MethPat.rowData works for 1-tuples", {
  expect_error(MethPat(rowData = granges(mt1)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'MTuples' object."))
})

test_that(".valid.MethPat.rowData works for 2-tuples", {
  expect_error(MethPat(rowData = granges(mt2)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'MTuples' object."))
})

test_that(".valid.MethPat.rowData works for 3-tuples", {
  expect_error(MethPat(rowData = granges(mt3)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'MTuples' object."))
})

test_that(".valid.MethPat.assays works for 1-tuples", {
  expect_error(MethPat(assays = list(), rowData = mt1), 
               paste0("Assay names must include all of: M, U"))
  # Extra assays are allowed
  ea <- c(assays(mp1), list('extraAssay' = 
                              matrix(1:20, ncol = 2, 
                                     dimnames = list(NULL, c('A', 'B')))))
  expect_is(MethPat(assays = ea, rowData = rowData(mp1)), "MethPat")
  # Assays must be non-negative (except extraAssays)
  a <- endoapply(ea, `-`, 10)
  expect_error(MethPat(assays = a, rowData = rowData(mp1)), 
               paste0("All counts of methylation patterns \\(stored in assays ", 
                      "slot\\) must be non-negative integers."))
  a <- ea
  a[['extraAssay']] <- a[['extraAssay']] - 100L
  expect_is(MethPat(assays = a, rowData = rowData(mp1)), "MethPat")
})

test_that(".valid.MethPat.assays works for 2-tuples", {
  expect_error(MethPat(assays = list(), rowData = mt2), 
               paste0("Assay names must include all of: MM, MU, UM, UU"))
  # Extra assays are allowed
  ea <- c(assays(mp2), list('extraAssay' = 
                              matrix(1:20, ncol = 2, 
                                     dimnames = list(NULL, c('A', 'B')))))
  expect_is(MethPat(assays = ea, rowData = rowData(mp2)), "MethPat")
  # Assays must be non-negative (except extraAssays)
  a <- endoapply(ea, `-`, 10)
  expect_error(MethPat(assays = a, rowData = rowData(mp2)), 
               paste0("All counts of methylation patterns \\(stored in assays ", 
                      "slot\\) must be non-negative integers."))
  a <- ea
  a[['extraAssay']] <- a[['extraAssay']] - 100L
  expect_is(MethPat(assays = a, rowData = rowData(mp2)), "MethPat")
})

test_that(".valid.MethPat.assays works for 3-tuples", {
  expect_error(MethPat(assays = list(), rowData = mt3), 
               paste0("Assay names must include all of: MMM, MMU, MUM, MUU, ", 
                      "UMM, UMU, UUM, UUU"))
  # Extra assays are allowed
  ea <- c(assays(mp3), list('extraAssay' = 
                              matrix(1:20, ncol = 2, 
                                     dimnames = list(NULL, c('A', 'B')))))
  expect_is(MethPat(assays = ea, rowData = rowData(mp3)), "MethPat")
  # Assays must be non-negative (except extraAssays)
  a <- endoapply(ea, `-`, 10)
  expect_error(MethPat(assays = a, rowData = rowData(mp3)), 
               paste0("All counts of methylation patterns \\(stored in assays ", 
                      "slot\\) must be non-negative integers."))
  a <- ea
  a[['extraAssay']] <- a[['extraAssay']] - 100L
  expect_is(MethPat(assays = a, rowData = rowData(mp3)), "MethPat")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("MethPat constructor")

test_that("MethPat constructor returns a valid MethPat object when m = 0", {
  expect_true(validObject(mp0))
})

test_that("MethPat constructor returns a valid object when m = 1", {
  expect_true(validObject(mp1))
})

test_that("MethPat constructor returns a valid object when m = 2", {
  expect_true(validObject(mp2))
})

test_that("MethPat constructor returns a valid object when m = 3", {
  expect_true(validObject(mp3))
})

test_that("MethPat constructor returns errors on bad input", {
  # TODO: None yet since the constructor doesn't check the input but relies on 
  # the validity methods.
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###
context("Combining MethPat objects")

test_that("cbind,MethPat-method works on good input", {
  expect_is(z <- cbind(mp1, mp1), "MethPat")
  expect_identical(dim(z), c(10L, 4L))
  expect_is(z <- cbind(mp1, mp1, mp1), "MethPat")
  expect_identical(dim(z), c(10L, 6L))
  expect_is(z <- cbind(mp2, mp2), "MethPat")
  expect_identical(dim(z), c(10L, 4L))
  expect_is(z <- cbind(mp2, mp2, mp2), "MethPat")
  expect_identical(dim(z), c(10L, 6L))
  expect_is(z <- cbind(mp3, mp3), "MethPat")
  expect_identical(dim(z), c(10L, 4L))
  expect_is(z <- cbind(mp3, mp3, mp3), "MethPat")
  expect_identical(dim(z), c(10L, 6L))
})

test_that("cbind,MethPat-method returns error on bad input", {
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp2), 
               "Cannot compare 'MTuples' objects of different 'size'.")
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp1[1:3]), 
               "'...' object ranges \\(rows\\) are not compatible")
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp1[10:1]), 
               "'...' object ranges \\(rows\\) are not compatible")
})

test_that("rbind,MethPat-method works on good input", {
  expect_is(z <- rbind(mp1, mp1), "MethPat")
  expect_identical(dim(z), c(20L, 2L))
  expect_is(z <- rbind(mp1, mp1, mp1), "MethPat")
  expect_identical(dim(z), c(30L, 2L))
  expect_is(z <- rbind(mp2, mp2), "MethPat")
  expect_identical(dim(z), c(20L, 2L))
  expect_is(z <- rbind(mp2, mp2, mp2), "MethPat")
  expect_identical(dim(z), c(30L, 2L))
  expect_is(z <- rbind(mp3, mp3), "MethPat")
  expect_identical(dim(z), c(20L, 2L))
  expect_is(z <- rbind(mp3, mp3, mp3), "MethPat")
  expect_identical(dim(z), c(30L, 2L))
})

test_that("rbind,MethPat-method returns error on bad input", {
  # TODO: Check error message is improved in new version of GenomicTuples
  expect_error(rbind(mp1, mp2), 
               "Cannot combine MTuples containing tuples of different 'size'.")
  mp1_ <- mp1
  colnames(mp1_) <- c('A', 'b')
  expect_error(rbind(mp1, mp1_), "'...' objects must have the same colnames")
})

test_that("combine,MethPat-method works for two MethPat objects", {
  # 1-tuples
  x <- mp1[1:2]
  y <- mp1[2:3]
  colnames(y) <- c('C', 'D')
  # Can't expect_identical because identical doesn't work on assay slot, 
  # because it is a refernece class (I think).
  expect_equal(combine(x, y), 
               MethPat(
                 assays = list(
                   M = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 
                                           2, 3)), ncol = 4, 
                              dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   U = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12, 13, 
                                           NA, 19, 18)), ncol = 4, 
                              dimnames = list(NULL, 
                                              c('A', 'B', 'C', 'D')))),
                 rowData = rowData(mp1)[1:3]
               )
  )
  # 2-tuples
  x <- mp2[1:2]
  y <- mp2[2:3]
  colnames(y) <- c('C', 'D')
  # Can't expect_identical because identical() returns false on assays slot, 
  # (I think this is because it is a reference class).
  expect_equal(combine(x, y), 
               MethPat(
                 assays = list(
                   MM = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 
                                            2, 3)), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MU = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12, 13, 
                                            NA, 19, 18)), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   UM = matrix(as.integer(c(30, 29, NA, 21, 22, NA, NA, 29, 28, 
                                            NA, 22, 23)), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UU = matrix(as.integer(c(40, 39, NA, 31, 32, NA, NA, 39, 38, 
                                            NA, 32, 33)), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D')))),
                 rowData = rowData(mp2)[1:3]
               )
  )
  # 3-tuples
  x <- mp3[1:2]
  y <- mp3[2:3]
  colnames(y) <- c('C', 'D')
  # Can't expect_identical because identical() returns false on assays slot, 
  # (I think this is because it is a reference class).
  expect_equal(combine(x, y), 
               MethPat(
                 assays = list(
                   MMM = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 
                                             2, 3)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MMU = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12, 
                                             13, NA, 19, 18)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   MUM = matrix(as.integer(c(30, 29, NA, 21, 22, NA, NA, 29, 
                                             28, NA, 22, 23)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MUU = matrix(as.integer(c(40, 39, NA, 31, 32, NA, NA, 39, 
                                             38, NA, 32, 33)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   UMM = matrix(as.integer(c(50, 49, NA, 41, 42, NA, NA, 49, 
                                             48, NA, 42, 43)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UMU = matrix(as.integer(c(60, 59, NA, 51, 52, NA, NA, 59, 
                                             58, NA, 52, 53)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UUM = matrix(as.integer(c(70, 69, NA, 61, 62, NA, NA, 69, 
                                             68, NA, 62, 63)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UUU = matrix(as.integer(c(80, 79, NA, 71, 72, NA, NA, 79, 
                                             78, NA, 72, 73)), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D')))),
                 rowData = rowData(mp3)[1:3]
               )
  )
})

test_that("combine,MethPat-method returns error on bad input", {
  x <- mp1[1:2]
  y <- mp1[2:3]
  expect_error(combine(x, y), 
               "Cannot combine 'MethPat' objects with duplicate colnames.")
  # TODO: Check error message is improved in new version of GenomicTuples
  expect_error(combine(mp1, mp2), 
               "Cannot combine MTuples containing tuples of different 'size'.")
  x <- mp1[1:2]
  y <- mp1[2:3]
  colnames(y) <- c('C', 'D')
  mcols(y) <- NULL
  # TODO: Write a more informative error message - might need to be specified 
  # for SummarizedExperiment.
  expect_error(combine(x, y), 
               "number of columns for arg 2 do not match those of first arg")
  x <- mp1[1:2]
  y <- mp1[2:3]
  colnames(y) <- c('C', 'D')
  assays(y) <- c(assays(y), list('extraAssay' = 
                                   matrix(1:4, ncol = 2, 
                                          dimnames = list(NULL, c('C', 'D')))))
  expect_error(combine(x, y), 
               "'MethPat' objects must all contain the same assays.")
  x <- mp3
  y <- mp3
  genome(y) <- 'mock2'
  expect_error(combine(x, y), 
               "sequences chr1, chr2, chr3 have incompatible genomes")
  y <- mp3
  colnames(y) <- c('C', 'D')
  seqlevelsStyle(y) <- 'NCBI'
  expect_warning(combine(x, y), 
               "The 2 combined objects have no sequence levels in common")
  y <- renameSeqlevels(y, c('chr1', '2', '3'))
  expect_identical(seqlevels(combine(x, y)), 
                   c('chr1', 'chr2', 'chr3', '2', '3'))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
context("MethPat getters")

test_that("SummarizedExperiment inherited getters work", {
  expect_identical(nrow(mp0), 0L)
  expect_identical(ncol(mp0), 0L)
  expect_identical(nrow(mp1), 10L)
  expect_identical(ncol(mp1), 2L)
  expect_identical(nrow(mp2), 10L)
  expect_identical(ncol(mp2), 2L)
  expect_identical(nrow(mp3), 10L)
  expect_identical(ncol(mp3), 2L)
  expect_identical(seqnames(mp1), mp1@rowData@seqnames)
  expect_identical(ranges(mp2), mp2@rowData@ranges)
  expect_identical(strand(mp3), mp3@rowData@strand)
  expect_identical(mcols(mp3), mp3@rowData@elementMetadata)
  expect_identical(elementMetadata(mp3), mp3@rowData@elementMetadata)
  expect_identical(seqinfo(mp3), mp3@rowData@seqinfo)
  expect_identical(seqlevels(mp3), seqlevels(mp3@rowData@seqinfo))
  expect_identical(seqlengths(mp3), seqlengths(mp3@rowData@seqinfo))
  expect_identical(isCircular(mp3), isCircular(mp3@rowData@seqinfo))
  expect_identical(genome(mp3), genome(mp3@rowData@seqinfo))
  expect_identical(seqlevelsStyle(mp3), seqlevelsStyle(mp3@rowData@seqinfo))
  # TODO: Notifiy Bioc-Devel that granges,SummarizedExperiment-method should 
  # return granges(rowData(x)) rather than rowData(x) since rowData may not be 
  # a GRanges object (e.g. might be a GTuples object)
  # expect_identical(granges(mp3), granges(mp3@rowData))
  expect_error(granges(mp3), "Not yet implemented")
})

test_that("methinfo getters work", {
  expect_identical(methinfo(mp1), methinfo(rowData(mp1)))
  expect_identical(methtype(mp1), methtype(rowData(mp1)))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###
context("MethPat splitting")

test_that("inherited split works", {
  # Split by integer
  mp3_s <- split(mp3, 1:10)
  expect_identical(length(mp3_s), 10L)
  expect_is(mp3_s, "SimpleList")
  expect_true(all(sapply(mp3_s, is, class = "MethPat")))
  # Split by Rle
  expect_message(mp3_s <- split(mp3, seqnames(mp3)), 
                 paste0("Note: method with signature ", 
                        sQuote("SummarizedExperiment#ANY"), " chosen for ", 
                        "function ", sQuote("split"), "."))
  expect_identical(length(mp3_s), 3L)
  expect_is(mp3_s, "SimpleList")
  expect_true(all(sapply(mp3_s, is, class = "MethPat")))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

context("MethPat setters")

test_that("SummarizedExperiment inherited setters work", {
  mp3_ <- mp3
  # TODO: Why isn't there a seqnames<-,SummarizedExperiment-method?
  expect_error(seqnames(mp3_) <- rev(seqnames(mp3)), 
               paste0("unable to find an inherited method for function ", 
                      sQuote("seqnames<-"), " for signature ", 
                      sQuote("\"MethPat\"")))
  mp3_ <- mp3
  ranges(mp3_) <- rev(ranges(mp3))
  expect_identical(ranges(mp3_), rev(ranges(mp3)))
  mp3_ <- mp3
  strand(mp3_) <- rev(strand(mp3))
  expect_identical(strand(mp3_), rev(strand(mp3)))
  mp3_ <- mp3
  mcols(mp3_) <- DataFrame(score = rev(mcols(mp3)$score))
  expect_identical(mcols(mp3_), DataFrame(score = rev(mcols(mp3)$score)))
  mp3_ <- mp3
  seqinfo(mp3_) <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                           seqlengths = c(10000L, 20000L, 15000L), 
                           isCircular = c(NA, NA, NA), 
                           genome = c("mock1", "mock1", "mock1"))
  expect_identical(seqinfo(mp3_), Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                                          seqlengths = c(10000L, 20000L, 
                                                         15000L), 
                                          isCircular = c(NA, NA, NA), 
                                          genome = c("mock1", "mock1", 
                                                     "mock1")))
  mp3_ <- mp3
  seqlevels(mp3_) <- c('chrI', 'chrII', 'chrIII')
  expect_identical(seqlevels(mp3_), c('chrI', 'chrII', 'chrIII'))
  mp3_ <- mp3
  seqlengths(mp3_) <- c(10000L, 20000L, 15000L)
  expect_identical(seqlengths(mp3_), c('chr1' = 10000L, 'chr2' = 20000L, 
                                       'chr3' = 15000L))
  mp3_ <- mp3
  isCircular(mp3_) <- c('chr1' = TRUE, 'chr2' = FALSE, 'chr3' = FALSE)
  expect_identical(isCircular(mp3_), c('chr1' = TRUE, 'chr2' = FALSE, 
                                       'chr3' = FALSE))
  mp3_ <- mp3
  genome(mp3_) <- 'foo'
  expect_identical(genome(mp3_), c('chr1' = 'foo', 'chr2' = 'foo', 
                                   'chr3' = 'foo'))
})

test_that("methinfo setters work", {
  methinfo(mp1) <- MethInfo(c('CG', 'CHG'))
  expect_identical(methinfo(mp1), MethInfo(c('CG', 'CHG')))
  methtype(mp1) <- c("CG", "CHG", "CHH")
  expect_identical(methtype(mp1), c("CG", "CHG", "CHH"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

test_that("IPD,Methpat-method works", {
  expect_error(IPD(mp0), "Cannot compute IPD from an empty 'MTuples'.")
  expect_error(IPD(mp1), 
               "It does not make sense to compute IPD when 'size' = 1.")
  expect_identical(IPD(mp2), IPD(mt2))
  expect_identical(IPD(mp3), IPD(mt3))
})

test_that("size,MethPat-method works", {
  expect_identical(size(mp0), NA_integer_)
  expect_identical(size(mp1), 1L)
  expect_identical(size(mp2), 2L)
  expect_identical(size(mp3), 3L)
})

test_that("tuples,MethPat-method works", {
  expect_identical(tuples(mp0), tuples(gt0))
  expect_identical(tuples(mp1), tuples(mt1))
  expect_identical(tuples(mp2), tuples(mt2))
  expect_identical(tuples(mp3), tuples(mt3))
})

test_that("tuples<-,MethPat-method works", {
  tuples(mp1) <- matrix(101:110, ncol = 1)
  expect_identical(tuples(mp1), 
                   matrix(101:110, ncol = 1, dimnames = list(NULL, 'pos1')))
  tuples(mp2) <- matrix(c(101:110, 102:111), ncol = 2)
  expect_identical(tuples(mp2), 
                   matrix(c(101:110, 102:111), ncol = 2, 
                          dimnames = list(NULL, c('pos1', 'pos2'))))
  tuples(mp3) <- matrix(c(101:110, 102:111, 103:112), ncol = 3)
  expect_identical(tuples(mp3), 
                   matrix(c(101:110, 102:111, 103:112), ncol = 3, 
                          dimnames = list(NULL, c('pos1', 'pos2', 'pos3'))))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###
context("MethPat subsetting")

# No tests yet. 
# Subsetting behaviour is entirely inherited via SummarizedExperiment.
