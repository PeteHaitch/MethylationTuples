### =========================================================================
### readMethtuple: Read methtuple .tsv files and create a MethPat object
### -------------------------------------------------------------------------
###

# TODO: Allow parallel reading of files or at least determine why it's not 
# possible. ANSWER: Reading-in in parallel is not a problem, rather, the 
# problem is with parallel merging of the files. This is complicated and so I 
# defer this for now.
# TODO: Profile.
# TODO: Add timing output?
# TODO: Compute approximate memory usage.
# TODO: Do argument checks before reading in the file, e.g., check seqinfo is a 
# valid seqinfo, so that functions fails fast+early.
# TODO: Will report 2 warnings if seqinfo is 'incorrect/insufficient'. The 
# second might better be suppressed, "2: In .Seqinfo.mergexy(x, y) :
# The 2 combined objects have no sequence levels in common. (Use
# suppressWarnings() to suppress this warning.)""

#' Read \code{.tsv} output files from \code{methtuple} software.
#'
#' @description
#' Read the \code{.tsv} output files from \code{methtuple} and construct a single 
#' \code{\link{MethPat}} object.  \code{methtuple} 
#' (\url{www.github.com/PeteHaitch/methtuple}) is Python software to extract 
#' methylation patterns at methylation loci from \code{BAM} files. 
#' All files should contain the same size m-tuples, e.g., all 2-tuples. All 
#' files should be mapped against the same reference genome, which is supplied 
#' as the \code{seqinfo} argument. 
#'
#' @param files The \code{.tsv} files created by \code{comethylation}. These 
#' files may be compressed with gzip or bzip2 (not yet implemented). All files 
#' must contain the same sized m-tuples. Files will be decompressed to a 
#' temporary directory via a call to \code{\link[base]{tempdir}}.
#' @param sample_names The sample names of each file. Must be unique.
#' @param methinfo A \code{\link{MethInfo}} object containing information about 
#' the the methylation loci in the \code{files}. This should be the minimal 
#' \code{MethInfo} object necessary to describe all methylation loci in all 
#' \code{files}. For example, if reading in two files, one with \code{CG} and 
#' one with \code{CHG} methylation, then the \code{methinfo} should be 
#' \code{MethInfo(c('CG', 'CHG'))}. \strong{NB:} It is generally recommended 
#' that you only store samples with the same \code{methinfo} in each 
#' \code{MethPat} object.
#' @param seqinfo A \code{\link[GenomeInfoDb]{Seqinfo}} object containing 
#' information about the reference genome of the samples. If none is supplied 
#' (\code{NULL}) then a bare-bones \code{\link[GenomeInfoDb]{Seqinfo}} object 
#' will be created containing only the \code{seqnames} inferred from the 
#' \code{files}.
#' @param verbose A \code{logical(1)} indicating whether messages about 
#' the reading of the data (via \code{data.table::\link[data.table]{fread}}) 
#' and about data coercion during construction of the \code{\link{MethPat}} 
#' object should be printed. Regardless of whether \code{verbose} is 
#' \code{TRUE} or \code{FALSE}, \code{readMethtuple} reports its progress via 
#' calls to \code{\link[base]{message}}; these can be suppressed by wrapped the 
#' call to \code{readMethtuple} in \code{\link[base]{suppressMesssages}}.
#' @param bpparam A \code{\link[BiocParallel]{bpparam}} object specifying the 
#' parallelisation strategy, if any. See below for a discussion of 
#' parallelisation options available with \code{readMethtuple}.
#' 
#' @section Parallelisation:
#' Parallelisation of \code{readMethtuple} is partially supported. Files 
#' may be decompressed in parallel but not read-in in parallel. 
#' Parallelisation uses the \pkg{BiocParallel} package. By default this uses a
#' \code{\link[BiocParallel]{MulticoreParam}()} instance or the user's 
#' preferred back-end if they have used \code{\link[BiocParallel]{register}}. 
#' Please consult the \pkg{BiocParallel} documentation for details on 
#' registering a parallel backend and parallelisation support available on 
#' different operating systems.
#'
#' @seealso \code{\link{MethPat}}
#' @return A \code{\link{MethPat}} object
#' @examples
#' # Using example data supplied with the MethylationTuples package
#' # Each file contains data on 20,000 tuples from a whole-genome 
#' # bisulfite-sequencing experiment on a human frontal cortex sample
#' # (http://www.ncbi.nlm.nih.gov/sra/?term=SRR949193).
#' 
#' # 1-tuples
#' x <- readMethtuple(system.file("extdata", "SRR949193_20k.rmdup.CG.1.tsv.gz", 
#' package = "MethylationTuples"), sample_names = "SRR949193_20k", 
#' methinfo = MethInfo('CG'), seqinfo = Seqinfo(seqnames = 'chr1', 
#' seqlengths = 249250621, isCircular = FALSE, genome = 'hg19'))
#' 
#' # 3-tuples
#' y <- readMethtuple(system.file("extdata", "SRR949193_20k.rmdup.CG.3.tsv.gz", 
#' package = "MethylationTuples"), sample_names = "SRR949193_20k", 
#' methinfo = MethInfo('CG'), seqinfo = Seqinfo(seqnames = 'chr1', 
#' seqlengths = 249250621, isCircular = FALSE, genome = 'hg19'))
#' 
#' # Can't mix 1-tuples and 3-tuples
#' \dontrun{
#' readMethtuple(c(system.file("extdata", "SRR949193_20k.rmdup.CG.1.tsv.gz", 
#' package = "MethylationTuples"), system.file("extdata", 
#' "SRR949193_20k.rmdup.CG.3.tsv.gz", package = "MethylationTuples")))
#' }
#' 
#' @export
readMethtuple <- function(files, 
                           sample_names = paste0('sample_', seq_along(files)), 
                           methinfo = MethInfo(), seqinfo = NULL, 
                           verbose = getOption('verbose'), bpparam = bpparam()) {
  
  # Check that there is a unique sample name for each filename
  if (length(sample_names) != length(files) | 
        length(sample_names) != length(unique(sample_names))) {
    stop("'sample_names' must be unique.")
  }
  
  # Check that the sample names don't contain the delimeter used in appending 
  # sample names to variables. While this could be a function argument, for now 
  # I simply choose a delimeter that is very unlikely to appear in sample names 
  # and include a check; this keeps a simpler interface to the function.
  DELIMETER <- '@@@'
  DELIMETER_REGEXP <- '@@@'
  if (any(grepl(pattern = DELIMETER_REGEXP, sample_names))) {
    stop("'sample_names' must not contain '", DELIMETER, "'.")
  }

  # TODO: Use R.utils::compressedFile() and associated methods.
  # Decompress files (if required).
  my_unzip <- function(files, verbose) {
    bplapply(files, function(file, verbose) {
      if (isGzipped(file)) {
        file <- gunzip(file, temporary = TRUE, remove = FALSE, skip = TRUE, 
                       overwrite = TRUE)
      } else if (isBzipped(file)) {
        file <- bunzip2(file, temporary= TRUE, remove = FALSE, skip = TRUE, 
                        overwrite = TRUE)
      }
      return(file)
    })
  }
  if (any(sapply(files, isGzipped)) || any(sapply(files, isBzipped))) {
    message("Decompressing files ...")
  }
  files <- unlist(my_unzip(files, verbose))
  
  # Read in file(s) serially and, if more than one file, merge these files.
  if (length(files) == 1L) {
    message(paste0("Reading file ", files))
  } else {
    message("Reading and merging files ...")
    message(paste0("\tReading ", files[1]))
  }
  mtsv <- fread(input = files[1], header = TRUE, verbose = verbose)
  keys <- c("chr", "strand", 
            grep(pattern = '^pos', x = names(mtsv), value = TRUE))
  setkeyv(x = mtsv, cols = keys, verbose = verbose)
  # Append the sample name to the names of the counts, e.g. 'M' -> 'M.sample_1'
  wn <- which(!(names(mtsv) %in% keys))
  setnames(mtsv, c(names(mtsv)[seq_len(ncol(mtsv))[-wn]], 
                   paste(names(mtsv)[wn], sample_names[1], sep = DELIMETER)))
  m <- sum(grepl(pattern = '^pos', x = names(mtsv)))
  strand <- unique(mtsv[, strand])
  if ('*' %in% strand && ('+' %in% strand || '-' %in% strand)) {
    warning(paste0("Some tuples are stranded ('strand' = '+' or '-') and ", 
                   "some are unstranded ('strand' = '*').\nIt is ", 
                   "recommended that all 'files' are identically stranded ", 
                   "or unstranded."))
  }
  
  if (length(files) > 1) {
    for(i in seq_along(files)[-1]) {
      message(paste0("\tReading ", files[i]))
      tmp <- fread(input = files[i], header = TRUE, verbose = verbose)
      strand <- unique(c(strand, unique(tmp[, strand])))
      if ('*' %in% strand && ('+' %in% strand || '-' %in% strand)) {
        warning(paste0("Some tuples are stranded ('strand' = '+' or '-') and ", 
                       "some are unstranded ('strand' = '*').\nIt is ", 
                       "recommended that all 'files' are identically stranded ", 
                       "or unstranded."))
      }
      message(paste0("\t\tMerging ", files[i]))
      setkeyv(x = tmp, cols = keys, verbose = verbose)
      mm <- sum(grepl(pattern = '^pos', x = names(tmp)))
      if (m != mm) {
        stop(paste0("'files' contain different sized tuples: ", m, " and ", mm))
      }
      setnames(tmp, c(names(tmp)[seq_len(ncol(tmp))[-wn]], 
                      paste(names(tmp)[wn], sample_names[i], sep = DELIMETER)))
      mtsv <- merge(mtsv, tmp, by = keys, all = TRUE)
    }
  }
  
  # Create the MethPat object
  message("Creating MethPat object ...")
  
  # Combine the counts of methylation patterns into the assays list
  assay_names <- unique(sapply(strsplit(grep(pattern = DELIMETER_REGEXP, 
                                             x = names(mtsv), value = TRUE),
                                        split = DELIMETER_REGEXP), '[[', 1))
  assays <- lapply(assay_names, function(an, mtsv) {
    pat <- paste0('^', an, '.')
    i <- grep(pattern = pat, x = names(mtsv), value = TRUE)
    sn <- sapply(strsplit(x = i, split = DELIMETER_REGEXP), '[[', 2)
    as.matrix(setnames(mtsv[, i, with = FALSE], sn))
  }, mtsv = mtsv)
  names(assays) <- assay_names
  
  sn <- unique(mtsv[, chr])
  if (is.null(seqinfo) || any(is.na(seqinfo@seqnames)) || 
        any(is.na(seqinfo@seqlengths)) || any(is.na(seqinfo@genome))) {
    warning(paste0("It is recommended that you supply a complete 'seqinfo' ", 
                   "including 'seqnames', 'seqlengths', 'isCircular' and ", 
                   "'genome'."))
    # Create a bare-bones Seqinfo object (can only infer seqnames from 'files').
    seqinfo <- Seqinfo(seqnames = sn)
  }
  sn_match <- match(sn, seqnames(seqinfo))
  if (any(is.na(sn_match))) {
    warning(paste0("'files' contain seqnames not found in 'seqinfo': ", 
                sn[which(is.na(sn_match))], "\n", 
                "Will try to create new seqinfo with all seqnames."))
    seqinfo <- merge(seqinfo, 
                     Seqinfo(seqnames = 
                               sn[which(is.na(sn_match))]))
  }
  
  rowRanges <- MTuples(
    gtuples = GTuples(seqnames = mtsv[, chr], 
                      strand = mtsv[, strand], 
                      tuples = as.matrix(mtsv[, grep(pattern = '^pos', 
                                                     x = colnames(mtsv)), 
                                              with = FALSE]), 
                      seqinfo = seqinfo), 
    methinfo = methinfo)
  
  new("MethPat", SummarizedExperiment(assays = assays, rowRanges = rowRanges,  
                                      verbose = verbose))
}
