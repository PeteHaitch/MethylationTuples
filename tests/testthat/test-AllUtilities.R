### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tests of utility functions
###
context("Utility functions")

test_that(".make_methpat_names works", {
  expect_identical(.make_methpat_names(1), c('M', 'U'))
  expect_identical(.make_methpat_names(2), c('MM', 'MU', 'UM', 'UU'))
  expect_identical(.make_methpat_names(3), c('MMM', 'MMU', 'MUM', 'MUU', 
                                             'UMM', 'UMU', 'UUM', 'UUU'))
})

test_that(".valid_methtype works", {
  
  expect_true(.valid_methtype('CG'))
  expect_true(.valid_methtype('CHG'))
  expect_true(.valid_methtype('CHH'))
  expect_true(.valid_methtype('CNN'))
  expect_true(.valid_methtype(c('CG', 'CHG')))
  expect_true(.valid_methtype(c('CG', 'CHH')))
  expect_true(.valid_methtype(c('CG', 'CNN')))
  expect_true(.valid_methtype(c('CHG', 'CHH')))
  expect_true(.valid_methtype(c('CHG', 'CNN')))
  expect_true(.valid_methtype(c('CHH', 'CNN')))
  expect_true(.valid_methtype(c('CG', 'CHG', 'CHH')))
  expect_true(.valid_methtype(c('CG', 'CHG', 'CNN')))
  expect_true(.valid_methtype(c('CG', 'CHH', 'CNN')))
  expect_true(.valid_methtype(c('CHG', 'CHH', 'CNN')))
  expect_true(.valid_methtype(c('CG', 'CHG', 'CHH', 'CNN')))  
  expect_false(.valid_methtype('CpG'))
  expect_false(.valid_methtype('CG/CHG'))
})

test_that(".stranded works", {
  expect_true(.stranded(GRanges('1', IRanges(1:3, 1:3), strand = "+")))
  expect_true(.stranded(GRanges('1', IRanges(1:3, 1:3), strand = "-")))
  expect_false(.stranded(GRanges('1', IRanges(1:3, 1:3), strand = "*")))
  expect_true(.stranded(GRanges('1', IRanges(1:3, 1:3), 
                                strand = c("+", "-", "+"))))
  expect_false(.stranded(GRanges('1', IRanges(1:3, 1:3), 
                                 strand = c("+", "*", "+"))))
  expect_false(.stranded(GRanges('1', IRanges(1:3, 1:3), 
                                 strand = c("-", "*", "-"))))
  expect_false(.stranded(GRanges('1', IRanges(1:3, 1:3), 
                                 strand = c("+", "*", "-"))))
})
