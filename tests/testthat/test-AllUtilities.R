### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tests of utility functions
###
context("Utility functions")

test_that(".makeMethPatNames works", {
  expect_identical(.makeMethPatNames(1), c("M", "U"))
  expect_identical(.makeMethPatNames(2), c("MM", "MU", "UM", "UU"))
  expect_identical(.makeMethPatNames(3), c("MMM", "MMU", "MUM", "MUU", 
                                           "UMM", "UMU", "UUM", "UUU"))
})

test_that(".validMethtype works", {
  
  expect_true(.validMethtype("CG"))
  expect_true(.validMethtype("CHG"))
  expect_true(.validMethtype("CHH"))
  expect_true(.validMethtype("CNN"))
  expect_true(.validMethtype(c("CG", "CHG")))
  expect_true(.validMethtype(c("CG", "CHH")))
  expect_true(.validMethtype(c("CG", "CNN")))
  expect_true(.validMethtype(c("CHG", "CHH")))
  expect_true(.validMethtype(c("CHG", "CNN")))
  expect_true(.validMethtype(c("CHH", "CNN")))
  expect_true(.validMethtype(c("CG", "CHG", "CHH")))
  expect_true(.validMethtype(c("CG", "CHG", "CNN")))
  expect_true(.validMethtype(c("CG", "CHH", "CNN")))
  expect_true(.validMethtype(c("CHG", "CHH", "CNN")))
  expect_true(.validMethtype(c("CG", "CHG", "CHH", "CNN")))  
  expect_false(.validMethtype("CpG"))
  expect_false(.validMethtype("CG/CHG"))
})

test_that(".isStranded works", {
  expect_true(.isStranded(GRanges("1", IRanges(1:3, 1:3), strand = "+")))
  expect_true(.isStranded(GRanges("1", IRanges(1:3, 1:3), strand = "-")))
  expect_false(.isStranded(GRanges("1", IRanges(1:3, 1:3), strand = "*")))
  expect_true(.isStranded(GRanges("1", IRanges(1:3, 1:3), 
                                strand = c("+", "-", "+"))))
  expect_false(.isStranded(GRanges("1", IRanges(1:3, 1:3), 
                                 strand = c("+", "*", "+"))))
  expect_false(.isStranded(GRanges("1", IRanges(1:3, 1:3), 
                                 strand = c("-", "*", "-"))))
  expect_false(.isStranded(GRanges("1", IRanges(1:3, 1:3), 
                                 strand = c("+", "*", "-"))))
})