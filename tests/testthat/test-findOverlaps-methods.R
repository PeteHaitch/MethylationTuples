# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps
###
context("findOverlaps,MethPat,MethPat-method")
mpq3 <- MethPat(assays = assays(mp3[1:9]), rowData = q3)
expect_identical(findOverlaps(mpq3, mpq3, type = 'any'), 
                 findOverlaps(q3, q3, type = 'any'))
expect_identical(findOverlaps(mpq3, mpq3, type = 'any', ignore.strand = TRUE), 
                 findOverlaps(q3, q3, type = 'any', ignore.strand = TRUE))
expect_identical(findOverlaps(mpq3, mpq3, type = 'start'), 
                 findOverlaps(q3, q3, type = 'start'))
expect_identical(findOverlaps(mpq3, mpq3, type = 'start', ignore.strand = TRUE), 
                 findOverlaps(q3, q3, type = 'start', ignore.strand = TRUE))
expect_identical(findOverlaps(mpq3, mpq3, type = 'end'), 
                 findOverlaps(q3, q3, type = 'end'))
expect_identical(findOverlaps(mpq3, mpq3, type = 'end', ignore.strand = TRUE), 
                 findOverlaps(q3, q3, type = 'end', ignore.strand = TRUE))
expect_identical(findOverlaps(mpq3, mpq3, type = 'within'), 
                 findOverlaps(q3, q3, type = 'within'))
expect_identical(findOverlaps(mpq3, mpq3, type = 'within', 
                              ignore.strand = TRUE), 
                 findOverlaps(q3, q3, type = 'within', ignore.strand = TRUE))
expect_identical(findOverlaps(mpq3, mpq3, type = 'equal'), 
                 findOverlaps(q3, q3, type = 'equal'))
expect_identical(findOverlaps(mpq3, mpq3, type = 'equal', ignore.strand = TRUE), 
                 findOverlaps(q3, q3, type = 'equal', ignore.strand = TRUE))