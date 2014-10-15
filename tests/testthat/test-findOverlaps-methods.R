# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps
###
context("findOverlaps,MethPat,MethPat-method")
mpmq3 <- MethPat(assays = assays(mp3[1:9]), rowData = mq3)
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'any'), 
                 findOverlaps(mq3, mq3, type = 'any'))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'any', ignore.strand = TRUE), 
                 findOverlaps(mq3, mq3, type = 'any', ignore.strand = TRUE))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'start'), 
                 findOverlaps(mq3, mq3, type = 'start'))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'start', ignore.strand = TRUE), 
                 findOverlaps(mq3, mq3, type = 'start', ignore.strand = TRUE))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'end'), 
                 findOverlaps(mq3, mq3, type = 'end'))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'end', ignore.strand = TRUE), 
                 findOverlaps(mq3, mq3, type = 'end', ignore.strand = TRUE))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'within'), 
                 findOverlaps(mq3, mq3, type = 'within'))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'within', 
                              ignore.strand = TRUE), 
                 findOverlaps(mq3, mq3, type = 'within', ignore.strand = TRUE))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'equal'), 
                 findOverlaps(mq3, mq3, type = 'equal'))
expect_identical(findOverlaps(mpmq3, mpmq3, type = 'equal', ignore.strand = TRUE), 
                 findOverlaps(mq3, mq3, type = 'equal', ignore.strand = TRUE))