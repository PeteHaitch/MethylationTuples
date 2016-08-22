# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("MTuples validity methods")

test_that("methinfo slot is checked", {
  expect_error(mt1@methinfo <- "kraken", 
               paste0("assignment of an object of class ", dQuote("character"), 
                      " is not valid for @", sQuote('methinfo'), " in an ", 
                      "object of class ", dQuote("MTuples"), "; is\\(value, ", 
                      "\"MethInfo\"\\) is not TRUE"))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("MTuples()")

test_that("MTuples() returns a valid MTuples object when size = 0", {
  expect_true(validObject(new("MTuples")))
  expect_true(validObject(MTuples()))
})

test_that("MTuples() returns a valid object when size > 0", {
  m <- 1:4
  lapply(m, function(mm) {
    expect_true(validObject(MTuples("chr1", matrix(seq.int(1, mm), ncol = mm))))
  })
})

test_that("MTuples() returns warnings on unexpected input", {
  expect_warning(MTuples("chr1", matrix(c(1.1, 2, 3), ncol = 1)), 
                 "Converting 'tuples' to integer mode")
})

test_that("MTuples() returns errors on bad input", {
  expect_error(MTuples("chr1", 1:10),
               "'tuples' must be an integer matrix")
  expect_error(MTuples("chr1", as.matrix(letters)), 
               "'tuples' must be an integer matrix")
  expect_error(MTuples("chr1", matrix(c(1, NA), ncol = 1)),
               "'NA' detected in 'tuples'")
  expect_error(MTuples("chr1", matrix(1L), methinfo = MethInfo('CpG')), 
               paste0("Invalid 'methtype'. Must be one or more of 'CG', 'CHG',", 
                      " 'CHH' or 'CN'"))
})

context("MTuplesFromGTuples()")

test_that("MTuplesFromGTuples() works when size = 0", {
  expect_true(validObject(MTuplesFromGTuples(gt0)))
})

test_that("MTuplesFromGTuples() works when size > 0", {
  expect_true(validObject(MTuplesFromGTuples(gt1)))
  expect_true(validObject(MTuplesFromGTuples(gt2)))
  expect_true(validObject(MTuplesFromGTuples(gt3)))
})

test_that("MTuplesFromGTuples() returns errors on bad input", {
  expect_error(MTuplesFromGTuples(as(gt2, "GRanges"), methinfo = MethInfo('CG')), 
               paste0("'gtuples' must be a 'GTuples' object"))
  expect_error(MTuplesFromGTuples(gt1, methinfo = MethInfo('CpG')), 
               paste0("Invalid 'methtype'. Must be one or more of 'CG', 'CHG',", 
                      " 'CHH' or 'CN'"))
})

context("MTuplesFromBSgenome()")

test_that(".OneTuplesFromDNAString() finds matches", {
  ir_cg <- .OneTuplesFromDNAString(dnastring, "CG")
  ir_chg <- .OneTuplesFromDNAString(dnastring, "CHG")
  ir_chh <- .OneTuplesFromDNAString(dnastring, "CHH")
  
  match_sequence <- mapply(function(ir, methtype) {
    n <- nchar(methtype[[1]])
    val <- logical(length(ir))
    for (i in seq_len(length(ir))) {
      if (as.logical(mcols(ir)$strand[i] == "+")) {
        val[i] <- dnastring[seq(start(ir)[i], 
                                start(ir)[i] + n - 1)] %in% 
          DNAStringSet(methtype)
      } else if (as.logical(mcols(ir)$strand[i] == "-")) {
        val[i] <- dnastring[seq(start(ir)[i] - n + 1, 
                                start(ir)[i])] %in% 
          reverseComplement(DNAStringSet(methtype))
      }
    }
    all(val)
  }, ir = list(ir_cg, ir_chg, ir_chh), 
  methtype = list(CG = c("CG"),
                  CHG = c("CAG", "CCG", "CTG", "CHG"),
                  CHH = c("CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", 
                          "CTC", "CTT", "CHH")))
  
  expect_true(all(match_sequence))
})

test_that(".OneTuplesFromDNAString() works on compound methtype", {
  
  ir_cg <- .OneTuplesFromDNAString(dnastring, "CG")
  ir_chg <- .OneTuplesFromDNAString(dnastring, "CHG")
  ir_chh <- .OneTuplesFromDNAString(dnastring, "CHH")
  ir_cg_chg <- .OneTuplesFromDNAString(dnastring, c("CG", "CHG"))
  ir_cg_chh <- .OneTuplesFromDNAString(dnastring, c("CG", "CHH"))
  ir_chg_chh <- .OneTuplesFromDNAString(dnastring, c("CHG", "CHH"))
  ir_cg_chg_chh <- .OneTuplesFromDNAString(dnastring, c("CG", "CHG", "CHH"))
  
  
  expect_identical(sort(.OneTuplesFromDNAString(dnastring, c("CG", "CHG"))),
                   sort(c(ir_cg, ir_chg)))
  expect_identical(sort(.OneTuplesFromDNAString(dnastring, c("CG", "CHH"))),
                   sort(c(ir_cg, ir_chh)))
  expect_identical(sort(.OneTuplesFromDNAString(dnastring, c("CHG", "CHH"))),
                   sort(c(ir_chg, ir_chh)))
  expect_identical(
    sort(.OneTuplesFromDNAString(dnastring, c("CG", "CHG", "CHH"))),
    sort(c(ir_cg, ir_chg, ir_chh)))
})

test_that(".MTuplesFromOneTuples() errors when size < 1", {
  expect_error(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                     size = 0,
                                     ignore.strand = FALSE), 
               "'size' must be a positive integer")
})
  

test_that(".MTuplesFromOneTuples() works when size = 1", {
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 1,
                                         ignore.strand = FALSE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 100L, 10L), 
                                             seq(6L, 105L, 10L)), 
                                           ncol = 1),
                           strand = strand(Rle(c("+", "-"), c(10, 10))),
                           seqinfo = seqinfo,
                           methinfo = mi))
})

test_that(".MTuplesFromOneTuples() works when size > 1", {
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 2,
                                         ignore.strand = FALSE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 90L, 10L), 
                                             seq(6L, 95L, 10L), 
                                             seq(11L, 100L, 10L),
                                             seq(16L, 105L, 10L)),
                                           ncol = 2),
                           strand = strand(Rle(c("+", "-"), c(9, 9))),
                           seqinfo = seqinfo,
                           methinfo = mi))
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 3,
                                         ignore.strand = FALSE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 80L, 10L), 
                                             seq(6L, 85L, 10L), 
                                             seq(11L, 90L, 10L),
                                             seq(16L, 95L, 10L),
                                             seq(21L, 100L, 10L),
                                             seq(26L, 105L, 10L)),
                                           ncol = 3),
                           strand = strand(Rle(c("+", "-"), c(8, 8))),
                           seqinfo = seqinfo,
                           methinfo = mi))
})

test_that(".MTuplesFromOneTuples() works when ignore.strand = TRUE", {
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 1,
                                         ignore.strand = TRUE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 100L, 5L)), 
                                           ncol = 1),
                           strand = strand(Rle("*", 20)),
                           seqinfo = seqinfo,
                           methinfo = mi))
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 2,
                                         ignore.strand = TRUE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 95L, 5L),
                                             seq(6L, 100L, 5L)), 
                                           ncol = 2),
                           strand = strand(Rle("*", 19)),
                           seqinfo = seqinfo,
                           methinfo = mi))
  expect_identical(.MTuplesFromOneTuples(one_tuples = one_tuples, 
                                         size = 3,
                                         ignore.strand = TRUE),
                   MTuples(seqnames = "chr1", 
                           tuples = matrix(c(seq(1L, 90L, 5L),
                                             seq(6L, 95L, 5L),
                                             seq(11L, 100L, 5L)), 
                                           ncol = 3),
                           strand = strand(Rle("*", 18)),
                           seqinfo = seqinfo,
                           methinfo = mi))
})

test_that(".MTuplesFromOneTuples() errors when size is too big", {
  expect_error(.MTuplesFromOneTuples(one_tuples = one_tuples,
                                     size = 100,
                                     ignore.strand = FALSE),
               "Cannot create tuples of size = 100 for seqnames = 'chr1")
})

test_that(".MTuplesFromOneTuples() errors on bad input", {
  expect_error(.MTuplesFromOneTuples(one_tuples = as(mt1, "GTuples"), 
                                     size = 2, 
                                     ignore.strand = TRUE),
               paste0("'one_tuples' must be an 'MTuples' object with 'size\\(", 
                      "one_tuples\\)' = 1"))
  expect_error(.MTuplesFromOneTuples(one_tuples = mt2, 
                                     size = 2, 
                                     ignore.strand = TRUE),
               paste0("'one_tuples' must be an 'MTuples' object with 'size\\(", 
                      "one_tuples\\)' = 1"))
})
  
test_that("MTuplesFromBSgenome() errors when size = 0", {
  expect_error(MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19,
                                   size = 0, 
                                   methinfo = mi,
                                   exclude = c(paste0("chr", c(1:9, "X", "Y")))),
               "'size' must be a positive integer")
})
test_that("MTuplesFromBSgenome() works when size > 0", {
  lapply(1:4, function(size) {
    x <- MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19,
                             size = size, 
                             methinfo = mi,
                             exclude = c(paste0("chr", c(1:9, "X", "Y")), "_"))
    expect_true(is(x, "MTuples"))
    expect_identical(size(x), size)
  })

})

test_that("MTuplesFromBSgenome() returns errors on bad input", {
  expect_error(
    MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19[["chrM"]],
                        size = 1,
                        methinfo = mi),
    "'bsgenome' must be a 'BSgenome' object")
  expect_error(MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19, 
                                   methinfo = "CG",
                                   size = 1),
               "'methinfo' must be a 'MethInfo' object")
  expect_error(MTuplesFromBSgenome(bsgenome = BSgenome.Hsapiens.UCSC.hg19, 
                                   methinfo = MethInfo("CN"),
                                   size = 1), 
               "'CN' methtype not yet supported")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###
context("MTuples coercion")

test_that("Coerction to GTuples works", {
  expect_is(as(mt0, "GTuples"), "GTuples")
  expect_is(as(mt1, "GTuples"), "GTuples")
  expect_is(as(mt2, "GTuples"), "GTuples")
  expect_is(as(mt3, "GTuples"), "GTuples")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning
###
context("MTuples updating and cloning")

test_that("update,MTuples-method works on all relevant slots", {
  mt1_update <- update(mt1, seqnames = rev(seqnames(mt1)))
  expect_identical(mt1_update, 
                   MTuples(seqnames = rev(seqnames(mt1)), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1), 
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1),
                           methinfo = methinfo(mt1)))
  mt1_update <- update(mt1, ranges = rev(ranges(mt1)))
  expect_identical(mt1_update, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = as.matrix(rev(tuples(mt1))), 
                           strand = strand(mt1), 
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  mt1_update <- update(mt1, strand = rev(strand(mt1)))
  expect_identical(mt1_update, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = rev(strand(mt1)),
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  mt1_update <- update(mt1, elementMetadata = DataFrame(score = Rle(0L, 10)))
  expect_identical(mt1_update, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1),
                           score = Rle(0L, 10), 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  seqinfo <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                     seqlengths = c(10000L, 20000L, 15000L), 
                     isCircular = c(NA, NA, NA), 
                     genome = c("mock1", "mock1", "mock1"))
  mt1_update <- update(mt1, seqinfo = seqinfo)
  expect_identical(mt1_update, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1),
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo,
                           methinfo = methinfo(mt1)))
  # metadata(gt1) is not the same as setting the metadata in the GTuples() 
  # constructor. This (somewhat confusing) behaviour is inherited from GRanges()
  mt1_update <- update(mt1, metadata = list("foo" = "bar"))
  mt1_metadata <- mt1
  metadata(mt1_metadata) <- list("foo" = "bar")
  expect_identical(mt1_update, mt1_metadata)
  methinfo <- MethInfo("CN")
  mt1_update <- update(mt1, methinfo = methinfo)
  expect_identical(mt1_update, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1),
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo))
  mt3_update <- update(mt3, ranges = IRanges(start(mt3) + 10L, end(mt3) + 10L), 
                       internalPos = mt3@internalPos + 10L)
  expect_identical(mt3_update, 
                   MTuples(seqnames = seqnames(mt3), 
                           tuples = unname(tuples(mt3)) + 10L, 
                           strand = strand(mt3), 
                           score = mcols(mt3)$score, 
                           seqinfo = seqinfo(mt3), 
                           methinfo = methinfo(mt3)))
})

test_that("clone,MTuples-method works", {
  mt1_clone <- GenomicRanges:::clone(mt1, seqnames = rev(seqnames(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(seqnames = rev(seqnames(mt1)), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1), 
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, ranges = rev(ranges(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = as.matrix(rev(tuples(mt1))), 
                           strand = strand(mt1), 
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, strand = rev(strand(mt1)))
  expect_identical(mt1_clone, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = rev(strand(mt1)), 
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  mt1_clone <- GenomicRanges:::clone(mt1, elementMetadata = 
                                       DataFrame(score = Rle(0L, 10)))
  expect_identical(mt1_clone, 
                   MTuples(seqnames = seqnames(mt1),
                           tuples = tuples(mt1), 
                           strand = strand(mt1),
                           score = Rle(0L, 10), 
                           seqinfo = seqinfo(mt1), 
                           methinfo = methinfo(mt1)))
  seqinfo <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                     seqlengths = c(10000L, 20000L, 15000L), 
                     isCircular = c(NA, NA, NA), 
                     genome = c("mock1", "mock1", "mock1"))
  mt1_clone <- GenomicRanges:::clone(mt1, seqinfo = seqinfo)
  expect_identical(mt1_clone, 
                   MTuples(seqnames = seqnames(mt1), 
                           tuples = tuples(mt1), 
                           strand = strand(mt1),
                           score = mcols(mt1)$score, 
                           seqinfo = seqinfo, 
                           methinfo = methinfo(mt1)))
  # metadata(mt1) is not the same as setting the metadata in the GTuples() 
  # constructor. This (somewhat confusing) behaviour is inherited from GRanges()
  mt1_clone <- GenomicRanges:::clone(mt1, metadata = list("foo" = "bar"))
  mt1_metadata <- mt1
  metadata(mt1_metadata) <- list("foo" = "bar")
  expect_identical(mt1_clone, mt1_metadata)
  mt3_clone <- GenomicRanges:::clone(mt3, 
                                     ranges = IRanges(start(mt3) + 10L, 
                                                      end(mt3) + 10L), 
                                     internalPos = mt3@internalPos + 10L)
  expect_identical(mt3_clone, 
                   MTuples(seqnames = seqnames(mt3),
                           tuples = unname(tuples(mt3)) + 10L, 
                           strand(mt3), 
                           score = mcols(mt3)$score, 
                           seqinfo = seqinfo(mt3), 
                           methinfo = methinfo(mt3)))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining and collapsing
###
context("Combining and collapsings MTuples")

test_that("c,MTuples-method works", {
  expect_identical(c(mt1), mt1)
  expect_identical(c(mt1[1:5], mt1[6:10]), mt1)
  expect_identical(c(mt2[1:5], mt2[6:10]), mt2)
  expect_identical(c(mt3[1:5], mt3[6:10]), mt3)
  expect_identical(do.call(c, list(mt3, mt3)), c(mt3, mt3))
  expect_error(c(mt3, as(mt3, "GTuples")), 
               paste0("unable to find an inherited method for function ", 
                      sQuote('methinfo'), " for signature ", 
                      sQuote('"GTuples"')))
  expect_error(c(mt3, granges(mt3)), 
               paste0("unable to find an inherited method for function ", 
                      sQuote('size'), " for signature ", sQuote('"GRanges"')))
  expect_error(c(mt3, mt2), 
               paste0("Cannot combine MTuples containing tuples of ", 
                      "different 'size'"))
  mt1_ <- mt1
  methinfo(mt1_) <- MethInfo("CHG")
  expect_identical(c(mt1, mt1_),
                   MTuplesFromGTuples(c(gt1, gt1), MethInfo(c("CG", "CHG"))))
})

test_that("strandCollapse,MTuples-method works on good input", {
   mt_pos <- MTuples(seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
                     tuples = matrix(c(10, 15, 100, 110, 150,
                                       15, 27, 110, 150, 154,
                                       27, 33, 150, 154, 166),
                                     ncol = 3),
                     strand = "+",
                     methinfo = MethInfo("CG"))
   mt_neg <- MTuples(seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
                     tuples = matrix(c(11, 16, 101, 111, 151,
                                       16, 28, 111, 151, 155,
                                       28, 34, 151, 155, 167),
                                     ncol = 3),
                     strand = "-",
                     methinfo = MethInfo("CG"))
  expect_identical(strandCollapse(mt_neg), mt_pos)
})

test_that("strandCollapse,MTuples-method errors on bad input", {
  expect_error(strandCollapse(mt1), "Object contains unstranded tuples")
  methtype(mt1) <- "CHG"
  expect_error(strandCollapse(mt1), 
               paste0("strandCollapse\\(\\) only supports 'MTuples' objects ",
                      "with the 'CG' methtype"))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
context("MTuples getters")

test_that("MTuples-specific getters work", {
  expect_identical(methinfo(mt1), MethInfo("CG"))
  expect_identical(methtype(mt1), "CG")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###
context("Splitting MTuples")

test_that("inherited split() works", {
  ## by integer
  mt2_s <- split(mt2, 1:10)
  expect_identical(length(mt2_s), 10L)
  expect_is(mt2_s, "MTuplesList")
  ## by Rle
  mt2_s <- split(mt2, seqnames(mt2))
  expect_identical(length(mt2_s), 3L)
  expect_is(mt2_s, "MTuplesList")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###
context("MTuples setters")

test_that("MTuples-specific setters work", {
  methinfo(mt1) <- MethInfo("CN")
  expect_identical(methinfo(mt1), MethInfo("CN"))
  methtype(mt1) <- "CHG"
  expect_identical(methtype(mt1), "CHG")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

# None required since tested in GenomicTuples package

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

# None required since tested in GenomicTuples package

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

# TODO: Not sure how to test this
