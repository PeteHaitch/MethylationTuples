# NB: Additional test data are defined in tests/testthat/fake_data/

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MethPat objects used in tests
###

mi <- MethInfo("CG")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GTuples objects used in tests
###

seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), rep(FALSE, 3), 
                   rep("mock1", 3))
# Empty GTuples object
gt0 <- GTuples(seqinfo = seqinfo)
# 1-tuples
gt1 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10), ncol = 1),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)
# 2-tuples
gt2 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10, 2:11), ncol = 2),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)
# 3-tuples
gt3 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10, 2:11, 3:12), ncol = 3),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)

# Construct a set of 3-tuples with "interesting" overlaps
q3 <- c(GTuples('chr1', matrix(as.integer(c(1, 1, 1, 3, 3, 5, 7, 9, 7)), 
                               ncol = 3), strand = '*'),
        GTuples('chr1', matrix(as.integer(c(1, 1, 1, 3, 3, 5, 7, 9, 7)), 
                               ncol = 3), strand = '+'),
        GTuples('chr1', matrix(as.integer(c(1, 1, 1, 3, 3, 5, 7, 9, 7)), 
                               ncol = 3), strand = '-'))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MTuples objects used in tests
###
mt0 <- MTuples(gt0, mi)
mt1 <- MTuples(gt1, mi)
mt2 <- MTuples(gt2, mi)
mt3 <- MTuples(gt3, mi)
mq3 <- MTuples(q3, mi)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MTuplesList objects used in tests
###
mtl0 <- MTuplesList(mt0)
mtl1 <- MTuplesList(mt1)
mtl2 <- MTuplesList(mt2)
mtl3 <- MTuplesList(mt3)
mql3 <- MTuplesList(mq3)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MethPat objects used in tests
###
# Empty MethPat object
mp0 <- MethPat()
# 1-tuples
mp1 <- MethPat(assays = 
                 SimpleList('M' = matrix(c(10:1, 1:10), ncol = 2,
                                         dimnames = list(NULL, c('A', 'B'))), 
                            'U' = matrix(c(11:20, 20:11), ncol = 2,
                                         dimnames = list(NULL, c('A', 'B')))), 
               rowRanges = mt1)
# 2-tuples
mp2 <- MethPat(assays = 
                 SimpleList('MM' = matrix(c(10:1, 1:10), ncol = 2,
                                         dimnames = list(NULL, c('A', 'B'))), 
                            'MU' = matrix(c(11:20, 20:11), ncol = 2,
                                         dimnames = list(NULL, c('A', 'B'))), 
                            'UM' = matrix(c(30:21, 21:30), ncol = 2,
                                          dimnames = list(NULL, c('A', 'B'))), 
                            'UU' = matrix(c(40:31, 31:40), ncol = 2,
                                          dimnames = list(NULL, c('A', 'B')))), 
               rowRanges = mt2)
# 3-tuples
mp3 <- MethPat(assays = 
                 SimpleList('MMM' = matrix(c(10:1, 1:10), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'MMU' = matrix(c(11:20, 20:11), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'MUM' = matrix(c(30:21, 21:30), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'MUU' = matrix(c(40:31, 31:40), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'UMM' = matrix(c(50:41, 41:50), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'UMU' = matrix(c(60:51, 51:60), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'UUM' = matrix(c(70:61, 61:70), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B'))), 
                            'UUU' = matrix(c(80:71, 71:80), ncol = 2,
                                           dimnames = list(NULL, c('A', 'B')))), 
               rowRanges = mt3)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Objects used in test of makeAdjacentPairsCpp and makeAllPairsCpp
###

gt <- GTuples(seqnames = Rle(c('chr1', 'chr3', 'chr2'), c(6, 6, 1)), 
              tuples = matrix(as.integer(c(1, 20, 30, 6, 16, 26, 1, 6, 11, 7, 
                                           17, 27, 10)), ncol = 1), 
              strand = Rle(c('+', '-', '+', '-', '+'), c(3, 3, 3, 3, 1)), 
              seqinfo = Seqinfo(seqnames = c('chr1', 'chr2', 'chr3')))
mt <- rev(MTuples(gt, MethInfo('CG')))
methpat <- MethPat(mt, assays = 
                     list("M" = matrix(1:13, ncol = 1, 
                                       dimnames = list(NULL, "A")), 
                          "U" = matrix(1:13, ncol = 1, 
                                       dimnames = list(NULL, "A"))))
methpat_order <- order(methpat)
methpat_rd_sorted <- sort(rowRanges(methpat))
