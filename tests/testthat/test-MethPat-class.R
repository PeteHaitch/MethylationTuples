# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("MethPat validity methods")

test_that(".valid.MethPat.rowData works for empty GTuples", {
  expect_error(MethPat(rowData = granges(gt0)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'GTuples' object."))
})

test_that(".valid.MethPat.rowData works for 1-tuples", {
  expect_error(MethPat(rowData = granges(gt1)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'GTuples' object."))
})

test_that(".valid.MethPat.rowData works for 2-tuples", {
  expect_error(MethPat(rowData = granges(gt2)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'GTuples' object."))
})

test_that(".valid.MethPat.rowData works for 3-tuples", {
  expect_error(MethPat(rowData = granges(gt3)), 
               paste0("'rowData' slot of a 'MethPat' object must be a ", 
                      "'GTuples' object."))
})

test_that(".valid.MethPat.assays works for 1-tuples", {
  expect_error(MethPat(assays = list(), rowData = gt1), 
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
  expect_error(MethPat(assays = list(), rowData = gt2), 
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
  expect_error(MethPat(assays = list(), rowData = gt3), 
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
  expect_is(cbind(mp1, mp1), "MethPat")
  expect_is(cbind(mp2, mp2), "MethPat")
  expect_is(cbind(mp3, mp3), "MethPat")
})

test_that("cbind,MethPat-method returns error on bad input", {
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp2), 
               "Cannot compare 'GTuples' objects of different 'size'.")
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp1[1:3]), 
               "'...' object ranges \\(rows\\) are not compatible")
  # TODO: Write a more informative error message.
  expect_error(cbind(mp1, mp1[10:1]), 
               "'...' object ranges \\(rows\\) are not compatible")
})

test_that("rbind,MethPat-method works on good input", {
  expect_is(rbind(mp1, mp1), "MethPat")
  expect_is(rbind(mp2, mp2), "MethPat")
  expect_is(rbind(mp3, mp3), "MethPat")
})

test_that("rbind,MethPat-method returns error on bad input", {
  # TODO: Write a more informative error message.
  expect_error(rbind(mp1, mp2), 
               "Cannot combine GTuples containing tuples of different 'size'.")
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
                   M = matrix(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 2, 3), 
                              ncol = 4, 
                              dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   U = matrix(c(11, 12, NA, 20, 19, NA, NA, 12, 13, NA, 19, 
                                18), ncol = 4, 
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
                   MM = matrix(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 2, 3), 
                               ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MU = matrix(c(11, 12, NA, 20, 19, NA, NA, 12, 13, NA, 19, 
                                 18), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   UM = matrix(c(30, 29, NA, 21, 22, NA, NA, 29, 28, NA, 22, 
                                 23), ncol = 4, 
                               dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UU = matrix(c(40, 39, NA, 31, 32, NA, NA, 39, 38, NA, 32,
                                 33), ncol = 4, 
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
                   MMM = matrix(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA, 2, 3), 
                                ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MMU = matrix(c(11, 12, NA, 20, 19, NA, NA, 12, 13, NA, 19, 
                                  18), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   MUM = matrix(c(30, 29, NA, 21, 22, NA, NA, 29, 28, NA, 22, 
                                  23), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   MUU = matrix(c(40, 39, NA, 31, 32, NA, NA, 39, 38, NA, 32,
                                  33), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))),
                   UMM = matrix(c(50, 49, NA, 41, 42, NA, NA, 49, 48, NA, 42,
                                  43), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UMU = matrix(c(60, 59, NA, 51, 52, NA, NA, 59, 58, NA, 52,
                                  53), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UUM = matrix(c(70, 69, NA, 61, 62, NA, NA, 69, 68, NA, 62,
                                  63), ncol = 4, 
                                dimnames = list(NULL, c('A', 'B', 'C', 'D'))), 
                   UUU = matrix(c(80, 79, NA, 71, 72, NA, NA, 79, 38, NA, 72,
                                  73), ncol = 4, 
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
  # TODO: Write a more informative error message.
  expect_error(combine(mp1, mp2), 
               "Cannot combine GTuples containing tuples of different 'size'.")
  x <- mp1[1:2]
  y <- mp1[2:3]
  colnames(y) <- c('C', 'D')
  mcols(y) <- NULL
  # TODO: Write a more informative error message.
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
  # TODO: Perhaps granges,SummarizedExperiment-method should return 
  # granges(rowData(x)) rather than rowData(x)
  expect_identical(granges(mp3), granges(mp3@rowData))  
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

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

# TODO: IPD, IPD, size, tuples, tuples<-

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

# TODO



