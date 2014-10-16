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
               rowData = mt1)
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
               rowData = mt2)
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
               rowData = mt3)
