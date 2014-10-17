# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("MTuples validity methods")

test_that("methinfo slot is checked", {
  expect_error(mt1@methinfo <- 'kraken', 
               paste0("assignment of an object of class ", dQuote("character"), 
                      " is not valid for @", sQuote('methinfo'), " in an ", 
                      "object of class ", dQuote("MTuples"), "; is\\(value, ", 
                      "\"MethInfo\"\\) is not TRUE"))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("MTuples constructor")

test_that("MTuples constructor returns a valid MTuples object when m = 0", {
  expect_true(validObject(new("MTuples")))
  expect_true(validObject(mt0))
})

test_that("MTuples constructor returns a valid object when m = 1", {
  expect_true(validObject(mt1))
})

test_that("MTuples constructor returns a valid object when m = 2", {
  expect_true(validObject(mt2))
})

test_that("MTuples constructor returns a valid object when m >= 3", {
  expect_true(validObject(mt3))
})

test_that("MTuples constructor returns warnings on unexpected input", {
  expect_warning(MTuples(GTuples('chr1', tuples = matrix(c(1.1, 2, 3), 
                                                         ncol = 1)), 
                         mi), 
                 "Converting 'tuples' to integer mode")
})

test_that("MTuples constructor returns errors on bad input", {
  expect_error(MTuples(GTuples('chr1', tuples = 1:10), mi),
               "'tuples' must be an integer matrix")
  expect_error(MTuples(GTuples('chr1', tuples = as.matrix(letters)), mi), 
               "'tuples' must be an integer matrix")
  expect_error(MTuples(GTuples('chr1', tuples = matrix(c(1, NA), ncol = 1)), 
                       mi),
               "'NA' detected in 'tuples'")
  expect_error(MTuples(mt1, MethInfo('CpG')), 
               paste0("Invalid 'methtype'. Must be one or more of 'CG', 'CHG',", 
                      " 'CHH' or 'CNN'"))
})


# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###
context("MTuples coercion")

test_that("Coerction to GTuples works", {
  expect_is(as(mt0, "GTuples"), "GTuples")
  expect_is(as(mt1, "GTuples"), "GTuples")
  expect_is(as(mt2, "GTuples"), "GTuples")
  expect_is(as(mt3, "GTuples"), "GTuples")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning
###
context("MTuples updating and cloning")

test_that("update works on all relevant slots", {
  mt1_update <- update(mt1, seqnames = rev(seqnames(mt1)))
  expect_identical(mt1_update, 
                   MTuples(GTuples(rev(seqnames(mt1)), tuples(mt1), 
                                   strand(mt1), score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_update <- update(mt1, ranges = rev(ranges(mt1)))
  expect_identical(mt1_update, 
                   MTuples(GTuples(seqnames(mt1), 
                                   as.matrix(rev(tuples(mt1))), 
                                   strand(mt1), score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_update <- update(mt1, strand = rev(strand(mt1)))
  expect_identical(mt1_update, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), 
                                   rev(strand(mt1)),
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_update <- update(mt1, elementMetadata = DataFrame(score = Rle(0L, 10)))
  expect_identical(mt1_update, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), strand(mt1),
                                   score = Rle(0L, 10), 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  seqinfo <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                     seqlengths = c(10000L, 20000L, 15000L), 
                     isCircular = c(NA, NA, NA), 
                     genome = c("mock1", "mock1", "mock1"))
  mt1_update <- update(mt1, seqinfo = seqinfo)
  expect_identical(mt1_update, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), strand(mt1),
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo), methinfo(mt1)))
  # metadata(gt1) is not the same as setting the metadata in the GTuples() 
  # constructor. This (somewhat confusing) behaviour is inherited from GRanges()
  mt1_update <- update(mt1, metadata = list("foo" = "bar"))
  mt1_metadata <- mt1
  metadata(mt1_metadata) <- list("foo" = "bar")
  expect_identical(mt1_update, mt1_metadata)
  methinfo <- MethInfo('CNN')
  mt1_update <- update(mt1, methinfo = methinfo)
  expect_identical(mt1_update, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), strand(mt1),
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), 
                           methinfo = methinfo))
  mt3_update <- update(mt3, ranges = IRanges(start(mt3) + 10L, end(mt3) + 10L), 
                       internalPos = mt3@internalPos + 10L)
  expect_identical(mt3_update, 
                   MTuples(GTuples(seqnames(mt3), unname(tuples(mt3)) + 10L, 
                                   strand(mt3), score = mcols(mt3)$score, 
                                   seqinfo = seqinfo(mt3)), methinfo(mt3)))
})

test_that("clone works", {
  mt1_clone <- GenomicRanges:::clone(mt1, seqnames = rev(seqnames(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(GTuples(rev(seqnames(mt1)), tuples(mt1), 
                                   strand(mt1), score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, ranges = rev(ranges(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(GTuples(seqnames(mt1), 
                                   as.matrix(rev(tuples(mt1))), strand(mt1), 
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, strand = rev(strand(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), 
                                   rev(strand(mt1)), 
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, elementMetadata = 
                                       DataFrame(score = Rle(0L, 10)))
  expect_identical(mt1_clone, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), strand(mt1),
                                   score = Rle(0L, 10), 
                                   seqinfo = seqinfo(mt1)), methinfo(mt1)))
  seqinfo <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                     seqlengths = c(10000L, 20000L, 15000L), 
                     isCircular = c(NA, NA, NA), 
                     genome = c("mock1", "mock1", "mock1"))
  mt1_clone <- GenomicRanges:::clone(mt1, seqinfo = seqinfo)
  expect_identical(mt1_clone, 
                   MTuples(GTuples(seqnames(mt1), tuples(mt1), strand(mt1),
                                   score = mcols(mt1)$score, 
                                   seqinfo = seqinfo), methinfo(mt1)))
  # metadata(mt1) is not the same as setting the metadata in the GTuples() 
  # constructor. This (somewhat confusing) behaviour is inherited from GRanges()
  mt1_clone <- GenomicRanges:::clone(mt1, metadata = list("foo" = "bar"))
  mt1_metadata <- mt1
  metadata(mt1_metadata) <- list("foo" = "bar")
  expect_identical(mt1_clone, mt1_metadata)
  mt3_clone <- GenomicRanges:::clone(mt3, ranges = IRanges(start(mt3) + 10L, 
                                                           end(mt3) + 10L), 
                                     internalPos = mt3@internalPos + 10L)
  expect_identical(mt3_clone, 
                   MTuples(GTuples(seqnames(mt3), unname(tuples(mt3)) + 10L, 
                                   strand(mt3), score = mcols(mt3)$score, 
                                   seqinfo = seqinfo(mt3)), methinfo(mt3)))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###
context("Combining MTuples")

test_that(".unlist_list_of_MTuples works", {
  expect_identical(.unlist_list_of_MTuples(list(mt0)), mt0)
  expect_identical(.unlist_list_of_MTuples(list(mt1)), mt1)
  expect_identical(.unlist_list_of_MTuples(list(mt2)), mt2)
  expect_identical(.unlist_list_of_MTuples(list(mt3)), mt3)
})

test_that("c works", {
  expect_identical(c(mt1[1:5], mt1[6:10]), mt1)
  expect_identical(c(mt2[1:5], mt2[6:10]), mt2)
  expect_identical(c(mt3[1:5], mt3[6:10]), mt3)
  expect_error(c(mt3, as(mt3, "GTuples")), 
               paste0("unable to find an inherited method for function ", 
                      sQuote('methinfo'), " for signature ", 
                      sQuote('"GTuples"')))
  expect_error(c(mt3, granges(mt3)), 
               paste0("unable to find an inherited method for function ", 
                      sQuote('size'), " for signature ", sQuote('"GRanges"')))
  expect_error(c(mt3, mt2), 
               paste0("Cannot combine MTuples containing tuples of ", 
                      "different 'size'"))
  mt1_ <- mt1
  methinfo(mt1_) <- MethInfo('CHG')
  expect_identical(c(mt1, mt1_),
    MTuples(c(as(mt1, "GTuples"), as(mt1, "GTuples")), 
                           MethInfo(c('CG', 'CHG'))))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
context("MTuples getters")

test_that("MTuples-specific getters work", {
  expect_identical(methinfo(mt1), MethInfo('CG'))
  expect_identical(methtype(mt1), 'CG')
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###
context("MTuples splitting")

test_that("inherited split works", {
  ## by integer
  mt2_s <- split(mt2, 1:10)
  expect_identical(length(mt2_s), 10L)
  expect_is(mt2_s, "MTuplesList")
  ## by Rle
  mt2_s <- split(mt2, seqnames(mt2))
  expect_identical(length(mt2_s), 3L)
  expect_is(mt2_s, "MTuplesList")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###
context("MTuples setters")

test_that("MTuples-specific setters work", {
  methinfo(mt1) <- MethInfo('CNN')
  expect_identical(methinfo(mt1), MethInfo('CNN'))
  methtype(mt1) <- 'CHG'
  expect_identical(methtype(mt1), 'CHG')
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

# None required since tested in GenomicTuples package

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

# None required since tested in GenomicTuples package

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

# TODO: Not sure how to test this. However, these tests are necessary due 
# to the showGTuples function being rather fragile.