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
a1 <- array(c(10:1, 1:10, 11:20, 20:11), 
            dim = c(10, 2, 2), 
            dimnames = list(NULL, c("A", "B"), c("M", "U")))
d1 <- DSArray::DSArray(a1, MARGIN = 2L)
h1 <- HDF5Array::HDF5Array(a1)
mpa1 <- MethPat(assays = a1, rowTuples = mt1)
mpd1 <- MethPat(assays = d1, rowTuples = mt1)
mph1 <- MethPat(assays = h1, rowTuples = mt1)

# 2-tuples
a2 <- array(c(10:1, 1:10, 11:20, 20:11, 30:21, 21:30, 40:31, 31:40),
            dim = c(10, 2, 4),
            dimnames = list(NULL, c("A", "B"), c("MM", "MU", "UM", "UU")))
d2 <- DSArray::DSArray(a2, MARGIN = 2L)
h2 <- HDF5Array::HDF5Array(a2)
mpa2 <- MethPat(assays = a2, rowTuples = mt2)
mpd2 <- MethPat(assays = d2, rowTuples = mt2)
mph2 <- MethPat(assays = h2, rowTuples = mt2)

# 2-tuples
a3 <- array(c(10:1, 1:10, 11:20, 20:11, 30:21, 21:30, 40:31, 31:40, 50:41, 
              41:50, 60:51, 51:60, 70:61, 61:70, 80:71, 71:80), 
            dim = c(10, 2, 8), 
            dimnames = list(NULL, c("A", "B"), c("MMM", "MMU", "MUM", "MUU", 
                                                 "UMM", "UMU", "UUM", "UUU")))
d3 <- DSArray::DSArray(a3, MARGIN = 2L)
h3 <- HDF5Array::HDF5Array(a3)
mpa3 <- MethPat(assays = a3, rowTuples = mt3)
mpd3 <- MethPat(assays = d3, rowTuples = mt3)
mph3 <- MethPat(assays = h3, rowTuples = mt3)
