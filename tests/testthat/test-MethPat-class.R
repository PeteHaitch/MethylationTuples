# # NB: Several objects used in testing are defined in
# # tests/testthat/helper-make-test-data.R
#
#

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("MethPat validity methods")

test_that(".valid.MethPat.rowTuples() checks class of rowTuples", {
  lapply(list(gt0, gt1, gt2, gt3), function(gt) {
    expect_error(MethPat(rowTuples = gt), 
                 paste0("'rowTuples' of 'MethPat' object must be an ",
                        "'MTuples' object."))
  })
})

test_that(".valid.MethPat.assays() checks for 'counts' assay", {
  lapply(list(mt1, mt2, mt3),  function(mt) {
    expect_error(MethPat(assays = list(), rowTuples = mt),
                 "Assays must include an element named 'counts'")
  })
})

test_that(".valid.MethPat.assays() allows for additional assays", {
  lapply(c(lmp1, lmp2, lmp3), function(mp) {
    ea <- c(assays(mp),
            list("extraAssay" = matrix(1, ncol = ncol(mp), nrow = nrow(mp),
                                       dimnames = list(NULL, colnames(mp)))))
    expect_is(MethPat(assays = ea,
                      rowTuples = rowTuples(mp),
                      colData = colData(mp)),
              "MethPat")
    # Assays must be non-negative (except extraAssays)
    a <- endoapply(ea, `-`, max(assay(mp, "counts")))
    expect_error(MethPat(assays = a,
                         rowTuples = rowTuples(mp),
                         colData = colData(mp)),
                 paste0("All 'counts' must be non-negative integers"))
    a <- ea
    a[["extraAssay"]] <- a[["extraAssay"]] - 100L
    expect_is(MethPat(assays = a,
                      rowTuples = rowTuples(mp),
                      colData = colData(mp)),
              "MethPat")
  })
})

test_that(".valid.MethPat.assays() checks slicenames(counts)", {
  lapply(c(lmp1, lmp2, lmp3), function(mp) {
    a <- assay(mp)
    dimnames(a)[[3]] <- letters[seq_len(dim(a)[3])]
    expect_error(MethPat(assays = a, 
                         rowTuples = rowTuples(mp), 
                         colData = colData(mp)), 
                 paste0("'counts' slicenames must be: '", 
                        paste0(.makeMethPatNames(size(mp)), 
                               collapse = "', '"), "'"))
  })
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("MethPat constructor")

# TODO: This feels redundant since the constructor calls the validity method
#       anyway
test_that("MethPat() returns a valid object", {
  # TODO: Actually call MethPat() rather than relying on objects created in 
  #       helper-make-test-data.R
  lapply(c(mp0, lmp1, lmp2, lmp3), function(mp1) {
    expect_true(validObject(mp1))
  })
})

test_that("MethPat() errors on bad input", {
  expect_error(MethPat(c1a, rowRanges = gt1),
               "'...' must not include 'rowRanges' or 'rowData'")
  expect_error(MethPat(c1a, rowData = gt1),
               "'...' must not include 'rowRanges' or 'rowData'")
  expect_error(MethPat(rowRanges = gt1),
               "'...' must not include 'rowRanges' or 'rowData'")
  expect_error(MethPat(rowData = gt1),
               "'...' must not include 'rowRanges' or 'rowData'")
  expect_error(MethPat(assays = c1a), 
               "Must supply an 'MTuples' object as the 'rowTuples' argument")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###
context("Combining MethPat objects")

test_that("cbind,MethPat-method works on good input", {
  lapply(c(lmp1, lmp2, lmp3), function(mp) {
    expect_is(z <- cbind(mp, mp), "MethPat")
    expect_identical(dim(z), c(10L, 4L))
    expect_is(z <- cbind(mp, mp, mp), "MethPat")
    expect_identical(dim(z), c(10L, 6L))
  })
})

test_that("cbind,MethPat-method returns error on bad input", {
  mapply(function(mp1, mp2) {
    expect_error(cbind(mp1, mp2),
                 "Cannot pcompare 'MTuples' objects of different 'size'.")
    expect_error(cbind(mp1, mp1[1:3]),
                 "'...' object ranges \\(rows\\) are not compatible")
    expect_error(cbind(mp1, mp1[10:1]),
                 "'...' object ranges \\(rows\\) are not compatible")
  }, mp1 = lmp1, mp2 = lmp2)
})

test_that("rbind,MethPat-method works on good input", {
  lapply(c(lmp1, lmp2, lmp3), function(mp) {
    expect_is(z <- rbind(mp, mp), "MethPat")
    expect_identical(dim(z), c(20L, 2L))
    expect_is(z <- rbind(mp, mp, mp), "MethPat")
    expect_identical(dim(z), c(30L, 2L))
  })
})

test_that("rbind,MethPat-method returns error on bad input", {
  mapply(function(mp1, mp2) {
    expect_error(rbind(mp1, mp2),
                 "Cannot combine MTuples containing tuples of different 'size'.")
    mp1_ <- mp1
    colnames(mp1_) <- c("A", "b")
    expect_error(rbind(mp1, mp1_), "'...' objects must have the same colnames")
  }, mp1 = lmp1, mp2 = lmp2)
})

# test_that("combine,MethPat-method works for two MethPat objects", {
#   # 1-tuples
#   x <- mp1[1:2]
#   y <- mp1[2:3]
#   colnames(y) <- c("C", "D")
#   # Can"t expect_identical because identical doesn"t work on assay slot,
#   # because it is a refernece class (I think).
#   expect_equal(combine(x, y),
#                MethPat(
#                  assays = list(
#                    M = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA,
#                                            2, 3)), ncol = 4),
#                    U = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12, 13,
#                                            NA, 19, 18)), ncol = 4)),
#                  rowTuples = rowTuples(mp1)[1:3],
#                  colData = DataFrame(row.names = c("A", "B", "C", "D"))
#                )
#   )
#   # 2-tuples
#   x <- mp2[1:2]
#   y <- mp2[2:3]
#   colnames(y) <- c("C", "D")
#   # Can"t expect_identical because identical() returns false on assays slot,
#   # (I think this is because it is a reference class).
#   expect_equal(combine(x, y),
#                MethPat(
#                  assays = list(
#                    MM = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA,
#                                             2, 3)), ncol = 4),
#                    MU = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12, 13,
#                                             NA, 19, 18)), ncol = 4),
#                    UM = matrix(as.integer(c(30, 29, NA, 21, 22, NA, NA, 29, 28,
#                                             NA, 22, 23)), ncol = 4),
#                    UU = matrix(as.integer(c(40, 39, NA, 31, 32, NA, NA, 39, 38,
#                                             NA, 32, 33)), ncol = 4)),
#                  rowTuples = rowRanges(mp2)[1:3],
#                  colData = DataFrame(row.names = c("A", "B", "C", "D"))
#                )
#   )
#   # 3-tuples
#   x <- mp3[1:2]
#   y <- mp3[2:3]
#   colnames(y) <- c("C", "D")
#   # Can"t expect_identical because identical() returns false on assays slot,
#   # (I think this is because it is a reference class).
#   expect_equal(combine(x, y),
#                MethPat(
#                  assays = list(
#                    MMM = matrix(as.integer(c(10, 9, NA, 1, 2, NA, NA, 9, 8, NA,
#                                              2, 3)), ncol = 4),
#                    MMU = matrix(as.integer(c(11, 12, NA, 20, 19, NA, NA, 12,
#                                              13, NA, 19, 18)), ncol = 4),
#                    MUM = matrix(as.integer(c(30, 29, NA, 21, 22, NA, NA, 29,
#                                              28, NA, 22, 23)), ncol = 4),
#                    MUU = matrix(as.integer(c(40, 39, NA, 31, 32, NA, NA, 39,
#                                              38, NA, 32, 33)), ncol = 4),
#                    UMM = matrix(as.integer(c(50, 49, NA, 41, 42, NA, NA, 49,
#                                              48, NA, 42, 43)), ncol = 4),
#                    UMU = matrix(as.integer(c(60, 59, NA, 51, 52, NA, NA, 59,
#                                              58, NA, 52, 53)), ncol = 4),
#                    UUM = matrix(as.integer(c(70, 69, NA, 61, 62, NA, NA, 69,
#                                              68, NA, 62, 63)), ncol = 4),
#                    UUU = matrix(as.integer(c(80, 79, NA, 71, 72, NA, NA, 79,
#                                              78, NA, 72, 73)), ncol = 4)),
#                  rowTuples = rowRanges(mp3)[1:3],
#                  colData <- DataFrame(row.names = c("A", "B", "C", "D"))
#                )
#   )
# })

# test_that("combine,MethPat-method returns error on bad input", {
#   x <- mp1[1:2]
#   y <- mp1[2:3]
#   expect_error(combine(x, y),
#                "Cannot combine 'MethPat' objects with duplicate colnames.")
#   # TODO: Check error message is improved in new version of GenomicTuples
#   expect_error(combine(mp1, mp2),
#                "Cannot combine MTuples containing tuples of different 'size'.")
#   x <- mp1[1:2]
#   y <- mp1[2:3]
#   colnames(y) <- c("C", "D")
#   mcols(y) <- NULL
#   # TODO: Write a more informative error message - might need to be specified
#   # for SummarizedExperiment.
#   expect_error(combine(x, y),
#                "number of columns for arg 2 do not match those of first arg")
#   x <- mp1[1:2]
#   y <- mp1[2:3]
#   colnames(y) <- c("C", "D")
#   assays(y) <- c(assays(y), list("extraAssay" =
#                                    matrix(1:4, ncol = 2,
#                                           dimnames = list(NULL, c("C", "D")))))
#   expect_error(combine(x, y),
#                "'MethPat' objects must all contain the same assays.")
#   x <- mp3
#   y <- mp3
#   genome(y) <- "mock2"
#   expect_error(combine(x, y),
#                "sequences chr1, chr2, chr3 have incompatible genomes")
#   y <- mp3
#   colnames(y) <- c("C", "D")
#   seqlevelsStyle(y) <- "NCBI"
#   expect_warning(combine(x, y),
#                  "The 2 combined objects have no sequence levels in common")
#   y <- renameSeqlevels(y, c("chr1", "2", "3"))
#   expect_identical(seqlevels(combine(x, y)),
#                    c("chr1", "chr2", "chr3", "2", "3"))
# })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
context("MethPat getters")

test_that("dim() works", {
  expect_identical(dim(mp0), c(0L, 0L))
  lapply(c(lmp1, lmp2, lmp3), function(mp) {
    expect_identical(dim(mp), c(10L, 2L))
  })
})

test_that("methinfo() works", {
  lapply(lmp1, function(mp1) {
    expect_identical(methinfo(mp1), methinfo(rowRanges(mp1)))
    expect_identical(methtype(mp1), methtype(rowRanges(mp1)))
  })
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###
context("MethPat splitting")

test_that("inherited split works", {
  # Split by integer
  lapply(lmp3, function(mp3) {
    mp3_s <- S4Vectors::split(mp3, 1:10)
    expect_identical(length(mp3_s), 10L)
    expect_is(mp3_s, "SimpleList")
    expect_true(all(sapply(mp3_s, is, class = "MethPat")))
    # Split by Rle
    mp3_s <- S4Vectors::split(mp3, seqnames(mp3))
    expect_identical(length(mp3_s), 3L)
    expect_is(mp3_s, "SimpleList")
    expect_true(all(sapply(mp3_s, is, class = "MethPat")))
  })
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

context("MethPat setters")

test_that("methinfo<- works", {
  lapply(lmp1, function(mp1) {
    methinfo(mp1) <- MethInfo(c("CG", "CHG"))
    expect_identical(methinfo(mp1), MethInfo(c("CG", "CHG")))
    methtype(mp1) <- c("CG", "CHG", "CHH")
    expect_identical(methtype(mp1), c("CG", "CHG", "CHH"))
  })
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

# TODO: Test
test_that("IPD,Methpat-method works", {
  expect_error(IPD(mp0), "Cannot compute IPD from an empty 'MTuples'.")
  lapply(lmp1, function(mp1) {
    expect_error(IPD(mp1),
                 "It does not make sense to compute IPD when 'size' = 1.")
  })
  lapply(c(lmp2, lmp3), function(mp){
    expect_identical(IPD(mp), IPD(mp))
  })
})

test_that("size,MethPat-method works", {
  mapply(function(mp, size) {
    expect_identical(size(mp), size)
  }, mp = c(mp0, lmp1, lmp2, lmp3),
  size = c(NA_integer_, rep(1L, 3), rep(2L, 3), rep(3L, 3)))
})

test_that("tuples,MethPat-method works", {
  mapply(function(mp, mt) {
    expect_identical(tuples(mp), tuples(mt))
  }, mp = c(mp0, lmp1, lmp2, lmp3),
  mt = c(list(mt0), replicate(3, mt1), replicate(3, mt2), replicate(3, mt3)))
})

test_that("tuples<-,MethPat-method works", {
  lapply(lmp1, function(mp1) {
    tuples(mp1) <- matrix(101:110, ncol = 1)
    expect_identical(tuples(mp1),
                     matrix(101:110, ncol = 1, dimnames = list(NULL, "pos1")))
  })
  lapply(lmp2, function(mp2) {
    tuples(mp2) <- matrix(c(101:110, 102:111), ncol = 2)
    expect_identical(tuples(mp2),
                     matrix(c(101:110, 102:111), ncol = 2,
                            dimnames = list(NULL, c("pos1", "pos2"))))
  })
  lapply(lmp3, function(mp3) {
    tuples(mp3) <- matrix(c(101:110, 102:111, 103:112), ncol = 3)
    expect_identical(tuples(mp3),
                     matrix(c(101:110, 102:111, 103:112), ncol = 3,
                            dimnames = list(NULL, c("pos1", "pos2", "pos3"))))
  })
})
