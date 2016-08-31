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
mt0 <- MTuplesFromGTuples(gt0, mi)
mt1 <- MTuplesFromGTuples(gt1, mi)
mt2 <- MTuplesFromGTuples(gt2, mi)
mt3 <- MTuplesFromGTuples(gt3, mi)
mq3 <- MTuplesFromGTuples(q3, mi)

one_tuples <- MTuples(seqnames = "chr1", 
                      tuples = matrix(c(seq(1L, 100L, 10L), 
                                        seq(6L, 105L, 10L)), 
                                      ncol = 1),
                      strand = Rle(c("+", "-"), c(10, 10)),
                      seqinfo = seqinfo,
                      methinfo = mi)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MTuplesList objects used in tests
###
mtl0 <- MTuplesList(mt0)
mtl1 <- MTuplesList(mt1)
mtl2 <- MTuplesList(mt2)
mtl3 <- MTuplesList(mt3)
mql3 <- MTuplesList(mq3)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DNAString objects used in tests
###

# A string containing all trinucleotides beginning with C, as well as 
# special methtype strings CHG, CHH, and CN
# the different methylation types separated by poly-A
c_tri_nucleotides <-
  data.frame(tri = apply(expand.grid("C", Biostrings::DNA_BASES, 
                                     Biostrings::DNA_BASES,
                                     stringsAsFactors = FALSE), 
                         1, paste, collapse = ""),
             stringsAsFactors = FALSE)
c_tri_nucleotides$rc_tri <- 
  as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(c_tri_nucleotides$tri)))
special_methtypes <- data.frame(tri = c("CHG", "CHH", "CN"),
                                rc_tri = c("CDG", "DDG", "NG"))
methtype_df <- rbind(c_tri_nucleotides, special_methtypes)
methtype_string <- paste(methtype_df$tri,
                         collapse = "NNN")
dnastring <- c(Biostrings::DNAString(methtype_string), 
               Biostrings::reverseComplement(
                 Biostrings::DNAString(methtype_string)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MethPat objects used in tests
###

c1a <- array(c(10:1, 1:10, 11:20, 20:11), 
             dim = c(10, 2, 2), 
             dimnames = list(NULL, c("s1", "s2"), c("M", "U")))
c1d <- DSArray::DSArray(c1a, 2)
c1h <- HDF5Array::HDF5Array(c1a)

c2a <- array(c(10:1, 1:10, 11:20, 20:11, 30:21, 21:30, 31:40, 40:31),
             dim = c(10, 2, 4), 
             dimnames = list(NULL, 
                             c("s1", "s2"), 
                             c("MM", "MU", "UM", "UU")))
c2d <- DSArray::DSArray(c2a, 2)
c2h <- HDF5Array::HDF5Array(c2a)

c3a <- array(c(10:1, 1:10, 11:20, 20:11, 30:21, 21:30, 31:40, 40:31,
               50:41, 41:50, 51:60, 60:51, 70:61, 61:70, 71:80, 80:71),
             dim = c(10, 2, 8), 
             dimnames = list(NULL, 
                             c("s1", "s2"),
                             c("MMM", "MMU", "MUM", "MUU", 
                               "UMM", "UMU", "UUM", "UUU")))
c3d <- DSArray::DSArray(c3a, 2)
c3h <- HDF5Array::HDF5Array(c3a)

cd <- DataFrame(row.names = c("s1", "s2"))

# Empty MethPat object
mp0 <- MethPat()

# 1-tuples
mp1a <- MethPat(assays = c1a, 
                rowTuples = mt1, 
                colData = cd)
mp1d <- MethPat(assays = c1d, 
                rowTuples = mt1, 
                colData = cd)
mp1h <- MethPat(assays = c1h,
                rowTuples = mt1,
                colData = cd)
lmp1 <- list(mp1a = mp1a, mp1d = mp1d, mp1h = mp1h)

# 2-tuples
mp2a <- MethPat(assays = c2a, 
                rowTuples = mt2, 
                colData = cd)
mp2d <- MethPat(assays = c2d, 
                rowTuples = mt2, 
                colData = cd)
mp2h <- MethPat(assays = c2h,
                rowTuples = mt2,
                colData = cd)
lmp2 <- list(mp2a = mp2a, mp2d = mp2d, mp2h = mp2h)

# 3-tuples
mp3a <- MethPat(assays = c3a, 
                rowTuples = mt3, 
                colData = cd)
mp3d <- MethPat(assays = c3d, 
                rowTuples = mt3, 
                colData = cd)
mp3h <- MethPat(assays = c3h,
                rowTuples = mt3,
                colData = cd)
lmp3 <- list(mp3a = mp3a, mp3d = mp3d, mp3h = mp3h)

