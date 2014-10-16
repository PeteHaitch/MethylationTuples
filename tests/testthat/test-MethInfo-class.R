# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

context("MethInfo getters")

test_that("methtype,MethInfo-method works", {
  expect_identical(methtype(MethInfo(NA_character_)), NA_character_)
  expect_identical(methtype(MethInfo(c('CG'))), c('CG'))
  expect_identical(methtype(MethInfo(c('CHG'))), c('CHG'))
  expect_identical(methtype(MethInfo(c('CHH'))), c('CHH'))
  expect_identical(methtype(MethInfo(c('CNN'))), c('CNN'))
  expect_identical(methtype(MethInfo(c('CG', 'CHG'))), c('CG', 'CHG'))
  expect_identical(methtype(MethInfo(c('CG', 'CHH'))), c('CG', 'CHH'))
  expect_identical(methtype(MethInfo(c('CG', 'CNN'))), c('CG', 'CNN'))
  expect_identical(methtype(MethInfo(c('CG', 'CHG', 'CHH'))), 
                   c('CG', 'CHG', 'CHH'))
  expect_identical(methtype(MethInfo(c('CG', 'CHG', 'CNN'))), 
                   c('CG', 'CHG', 'CNN'))
  expect_identical(methtype(MethInfo(c('CG', 'CHH', 'CNN'))), 
                   c('CG', 'CHH', 'CNN'))
  expect_identical(methtype(MethInfo(c('CG', 'CHG', 'CHH', 'CNN'))), 
                   c('CG', 'CHG', 'CHH', 'CNN'))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

test_that("MethInfo validity works", {
  expect_error(MethInfo(NA), 
               "invalid object for slot \"methtype\" in class \"MethInfo\":")
  expect_error(MethInfo("CpG"), "Invalid 'methtype'")
  expect_error(MethInfo("CG/CHG"), "Invalid 'methtype'")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

test_that("MethInfo constructor works", {
  expect_is(MethInfo('CG'), 'MethInfo')
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###
test_that("methtype<-,MethInfo-method works", {
  methtype(mi) <- 'CHG'
  expect_identical(mi, MethInfo('CHG'))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show.
###

test_that("summary,MethInfo-method works", {
  expect_identical(summary(mi), "CG methylation type")
  expect_identical(summary(MethInfo()), "NA methylation type")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###
test_that("merge,MethInfo-method works", {
  expect_identical(merge(MethInfo('CHG'), MethInfo('CG')), 
                   MethInfo(c('CG', 'CHG')))
  expect_identical(merge(mi, mi), mi)
  expect_identical(merge(mi, MethInfo()), MethInfo())
})
