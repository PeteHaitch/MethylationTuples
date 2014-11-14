### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readMethtuple
###
context("readMethtuple")

test_that("Returns MethPat object on good input", {
  expect_is(suppressMessages(
    readMethtuple('fake_data/good_fake_1.tsv', 
                   seqinfo = seqinfo)), 
    "MethPat")
  expect_is(suppressMessages(
    readMethtuple('fake_data/good_fake_3.tsv', 
                   seqinfo = seqinfo)), 
    "MethPat")
})

test_that("Returns warning on bad input", {
  expect_warning(suppressMessages(
    readMethtuple('fake_data/good_fake_1.tsv')), 
    paste0("It is recommended that you supply a complete ", 
           "'seqinfo' including 'seqnames', 'seqlengths', ", 
           "'isCircular' and 'genome'."))
  expect_warning(suppressMessages(
    readMethtuple('fake_data/bad_fake_1.tsv', 
                   seqinfo = Seqinfo('chr1'))),
    paste0("It is recommended that you supply a complete ", 
           "'seqinfo' including 'seqnames', 'seqlengths', ", 
           "'isCircular' and 'genome'."))
  expect_warning(suppressMessages(
    readMethtuple('fake_data/bad_fake_1.tsv', 
                   seqinfo = seqinfo)), 
    paste0("Some tuples are stranded \\('strand' = '\\+' or '-'\\) and some ", 
           "are unstranded \\('strand' = '\\*'\\)"))
  expect_warning(suppressMessages(
    readMethtuple('fake_data/bad_fake_3.tsv', 
                   seqinfo = seqinfo)), 
    paste0("Some tuples are stranded \\('strand' = '\\+' or '-'\\) and some ", 
           "are unstranded \\('strand' = '\\*'\\)"))
})

test_that("sample_names are checked", {
  expect_error(
    suppressMessages(readMethtuple(
      c('fake_data/good_fake_1.tsv', 'fake_data/good_fake_1.tsv'),
      sample_names = c('a', 'a'), seqinfo = seqinfo)), 
    "'sample_names' must be unique.")
  expect_error(
    suppressMessages(readMethtuple(
      'fake_data/good_fake_1.tsv', sample_names = 'a@@@', 
      seqinfo = seqinfo)),
    "'sample_names' must not contain '@@@'.")
})

test_that("tuple sizes are checked", {
  expect_error(suppressMessages(readMethtuple(
    c('fake_data/good_fake_1.tsv', 'fake_data/good_fake_3.tsv'), 
    sample_names = c('a', 'b'), seqinfo = seqinfo)), 
    "'files' contain different sized tuples: 1 and 3."
  )
})

test_that("bare-bones seqinfo is constructed", {
  expect_identical(seqinfo(
    suppressMessages(suppressWarnings(
      readMethtuple('fake_data/good_fake_1.tsv')))), 
    Seqinfo('chr1'))
})

test_that("seqinfo is sanity checked", {
  expect_warning(
    suppressMessages(
      z <- readMethtuple('fake_data/good_fake_1.tsv', 
                          seqinfo = Seqinfo('chr2', 1000, FALSE, 'mock'))), 
    "'files' contain seqnames not found in 'seqinfo': chr1")
  expect_identical(seqinfo(z), 
                   suppressWarnings(merge(Seqinfo('chr2', 1000, FALSE, 'mock'), 
                                          Seqinfo('chr1'))))
})


# 
# files, 
# sample_names = paste0('sample_', seq_along(files)), 
# methinfo = MethInfo(), seqinfo = NULL, 
# verbose = getOption('verbose'), bpparam = bpparam()