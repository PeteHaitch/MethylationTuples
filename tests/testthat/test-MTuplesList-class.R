### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("MTuplesList constructor")

test_that("MTuplesList constructor returns a valid object with no arguments", {
  expect_true(validObject(MTuplesList()))
})

test_that("MTuplesList constructor works for different numbers of arguments", {
  expect_true(validObject(MTuplesList(mt1)))
  expect_true(validObject(MTuplesList(mt1, mt1)))
  expect_true(validObject(MTuplesList(mt1, mt1, mt1)))
})

test_that("MTuplesList constructor returns a valid object when m = 0", {
  expect_true(validObject(mtl0))
})

test_that("MTuplesList constructor returns a valid object when m = 1", {
  expect_true(validObject(mtl1))
})

test_that("MTuplesList constructor returns a valid object when m = 2", {
  expect_true(validObject(mtl2))
})

test_that("MTuplesList constructor returns a valid object when m >= 3", {
  expect_true(validObject(mtl3))
})

test_that("MTuplesList constructor returns errors on bad input", {
  expect_error(MTuplesList(mt1, mt2), 
               "all MTuples in '...' must have the same 'size'")
  gr <- granges(gt2)
  expect_error(MTuplesList(mt2, gt2), 
               "all elements in '...' must be MTuples objects")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters/setters
###
context("MTuplesList-specific accessors")

test_that("methinfo works", {
  expect_identical(methinfo(mtl1), mi)
})

test_that("methinfo<- works", {
  methinfo(mtl1) <- MethInfo(c("CHG", "CHH"))
  expect_identical(methinfo(mtl1), MethInfo(c("CHG", "CHH")))
})

test_that("methtype works", {
  expect_identical(methtype(mtl1), "CG")
})

test_that("methtype<- works", {
  methtype(mtl1) <- c("CNN")
  expect_identical(methtype(mtl1), "CNN")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###
context("MTuplesList coercion")

test_that("Coercion to GTuplesList works", {
  expect_is(as(mtl0, "GTuplesList"), "GTuplesList")
  expect_is(as(mtl1, "GTuplesList"), "GTuplesList")
  expect_is(as(mtl2, "GTuplesList"), "GTuplesList")
  expect_is(as(mtl3, "GTuplesList"), "GTuplesList")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deconstruction/reconstruction of a GTuplesList into/from a GTuples
### object.
###
### For internal use only (not exported).
###
context("Desconstruction/reconstruction of MTuplesList")
test_that("Returns identical", {
  expect_identical(GenomicRanges:::reconstructGRLfromGR(
    GenomicRanges:::deconstructGRLintoGR(mtl0), mtl0), mtl0)
  expect_identical(GenomicRanges:::reconstructGRLfromGR(
    GenomicRanges:::deconstructGRLintoGR(mtl1), mtl1), mtl1)
  expect_identical(GenomicRanges:::reconstructGRLfromGR(
    GenomicRanges:::deconstructGRLintoGR(mtl2), mtl2), mtl2)
  expect_identical(GenomicRanges:::reconstructGRLfromGR(
    GenomicRanges:::deconstructGRLintoGR(mtl3), mtl3), mtl3)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

# TODO: Not sure how to test this. However, these tests are necessary due 
# to the showGTuples function being rather fragile.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from MTuples to MTuplesList with extractList() and family.
###
context("GTuplesList relistToClass")

test_that("relistToClass works", {
  expect_identical(relistToClass(mt0), "MTuplesList")
  expect_identical(relistToClass(mt1), "MTuplesList")
  expect_identical(relistToClass(mtl0), "SimpleList")
  expect_identical(relistToClass(mtl1), "SimpleList")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### apply methods
###
context("MTuplesList apply methods")

test_that("MTuplesList endoapply works", {
  expect_identical(mtl0, endoapply(mtl0, function(x) {x}))
  expect_identical(mtl1, endoapply(mtl1, function(x) {x}))
  expect_identical(mtl2, endoapply(mtl2, function(x) {x}))
  expect_identical(mtl3, endoapply(mtl3, function(x) {x}))
})

