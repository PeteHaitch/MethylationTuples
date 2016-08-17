### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tests of utility functions
###
context("Utility functions")

test_that(".validMethtype works", {
  expect_true(.validMethtype("CG"))
  expect_true(.validMethtype("CHG"))
  expect_true(.validMethtype("CHH"))
  expect_true(.validMethtype("CN"))
  expect_true(.validMethtype(c("CG", "CHG")))
  expect_true(.validMethtype(c("CG", "CHH")))
  expect_true(.validMethtype(c("CG", "CN")))
  expect_true(.validMethtype(c("CHG", "CHH")))
  expect_true(.validMethtype(c("CHG", "CN")))
  expect_true(.validMethtype(c("CHH", "CN")))
  expect_true(.validMethtype(c("CG", "CHG", "CHH")))
  expect_true(.validMethtype(c("CG", "CHG", "CN")))
  expect_true(.validMethtype(c("CG", "CHH", "CN")))
  expect_true(.validMethtype(c("CHG", "CHH", "CN")))
  expect_true(.validMethtype(c("CG", "CHG", "CHH", "CN")))  
  expect_false(.validMethtype("CpG"))
  expect_false(.validMethtype("CG/CHG"))
  expect_false(.validMethtype(c(NA_character_, NA_character_)))
})
