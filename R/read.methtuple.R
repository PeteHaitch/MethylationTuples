### =========================================================================
### read.methtuple: Read methtuple .tsv files and create a MethPat object
### -------------------------------------------------------------------------
###

# TODO: Add bzip2 support. Do this by adding a bzip2 function to R.utils.
# TODO: Improved support of compressed files, e.g., better recognition of file 
# type and decompression to temporary folder.
# TODO: Allow parallel reading of files or at least determine why it's not 
# possible.
# TODO: Make verbose argument a verbosity argument with multiple levels. 
# 0 = no output, 1 = progress reports, 2 = full verbosity.
# TODO: Should I print/cat/message progress reports
# TODO: unit tests
# TODO: Use a different delimiter to '.' when adding sample names to count 
# names. Otherwise a sample name with '.' in it will break the code. Whatever 
# delimiter I end up using, add a check that this delimiter is not present in 
# any of the sample_names.
# TODO: Profile.
# TODO: Add timing output?
# TODO: Compute approximate memory usage.
# TODO: Check that all files only contain either strand %in% '*' or strand 
# %in% c('+', '-'); otherwise the resulting object is a bit of a mess.

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
#' must contain the same sized m-tuples.
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
#' information about the reference genome of the samples.
#' @param verbose A \code{\link{logical(1)}} indicating whether messages about 
#' the reading of the data (with \code{data.table::\link[data.table]{fread}}) 
#' and about data coercion during construction should be printed.
#' @param bpparam A \code{\link[BiocParallel]{bpparam}} object specifying the 
#' parallelisation strategy, if any. See below for a discussion of 
#' parallelisation options available with \code{read.methtuple}.
#' 
#' @section Parallelisation:
#' Parallelisation of \code{read.methtuple} is partially supported. Files 
#' may be decompressed in parallel but not read-in in parallel. 
#' Parallelisation uses the \pkg{BiocParallel} package. By default this uses a
#' \code{\link[BiocParallel]{MulticoreParam()}} instance or the user’s 
#' preferred back-end if they have used \code{\link[BiocParallel]{register}}. 
#' Please consult the \pkg{BiocParallel} documentation for details on 
#' registering a parallel backend and parallelisation support available on 
#' different operating systems.
#' 
#' @note The compression of \code{files} is determined by the file extensions: 
#' \code{.gz} for compressed with gzip and \code{.bz2} for compressed with 
#' bzip2. In all other cases the file is assumed to be uncompressed.
#' 
#' Files will be decompressed to a temporary directory.
#'
#' @seealso \code{\link{MethPat}}
#' @return A \code{\link{MethPat}} object
#' @examples
#' ## TODO
#' 
#' @export
read.methtuple <- function(files, 
                           sample_names = paste0('sample_', seq_along(files)), 
                           methinfo = MethInfo(), seqinfo = Seqinfo(), 
                           verbose = FALSE, bpparam = bpparam()) {
  
  # Check that there is a unique sample name for each filename
  if (length(sample_names) != length(files) | 
        length(sample_names) != length(unique(sample_names))) {
    stop("Each file must have a unique sample name.")
  }
  
  # Decompress files (if required).
  if (any(grepl("\\.gz$", files) || any(grepl("\\.bz2$", files)))) {
    print("Decompressing files...")
    
    my_unzip <- function(files, verbose) {
      bplapply(files, function(file, verbose) {
        if (grepl("\\.gz$", file)){
          file <- gunzip(file, temporary = TRUE, remove = FALSE)
        } else if (grepl("\\.bz2$", file)){
          stop("Sorry, can't yet handle ", sQuote('bzip2'), 
               " compressed files.")
          # TODO: Uncomment once I have a way to deal with bzip files
          # file <- bzip(file, temporary = TRUE, remove = FALSE)
        } else {
          # Nothing to do!
        }
        return(file)
      })
    }
    files <- unlist(my_unzip(files, verbose))
  } else {
    # Nothing to do!
  }
  
  # Read in file(s) serially and, if more than one file, merge these files.
  if (length(files) == 1L) {
    print(paste0("Reading file ", files))
  } else {
    print(paste0("Reading file ", files[1]))
  }
  mtsv <- fread(input = files[1], header = TRUE, verbose = verbose)
  keys <- c("chr", "strand", 
            grep(pattern = '^pos', x = names(mtsv), value = TRUE))
  setkeyv(x = mtsv, cols = keys, verbose = verbose)
  # Append the sample name to the names of the counts, e.g. 'M' -> 'M.sample_1'
  wn <- which(!(names(mtsv) %in% keys))
  setnames(mtsv, c(names(mtsv)[seq_len(ncol(mtsv))[-wn]], 
                   paste(names(mtsv)[wn], sample_names[1], sep = ".")))
  if (length(files) > 1) {
    for(i in seq_along(files)[-1]) {
      print(paste0("Reading and merging ", files[i]))
      ## TODO: try-catch the merge in case files with different sized m-tuples
      ## are supplied by the user because the error message supplied by 
      ## merge is a bit hard to understand.
      tmp <- fread(input = x$file[i], header = TRUE, verbose = verbose)
      setkeyv(x = tmp, cols = keys, verbose = verbose)
      wn <- which(!(names(tmp) %in% keys))
      setnames(tmp, c(names(tmp)[seq_len(ncol(tmp))[-wn]], 
                      paste(names(tmp)[wn], sample_names[i], sep = ".")))
      mtsv <- merge(mtsv, tmp, by = keys, all = TRUE)
    }
  }
  
  # Combine the counts into the assays list
  print("Creating assays...")
  assay_names <- unique(sapply(strsplit(grep(pattern = '\\.', x = names(mtsv), 
                                             value = TRUE),
                                        split = '\\.'), '[[', 1))
  assays <- lapply(assay_names, function(an, mtsv) {
    pat <- paste0('^', an, '.')
    i <- grep(pattern = pat, x = names(mtsv), value = TRUE)
    sn <- sapply(strsplit(x = i, split = '\\.'), '[[', 2)
    as.matrix(setnames(mtsv[, i, with = FALSE], sn))
  }, mtsv = mtsv)
  names(assays) <- assay_names
  
  rowData <- MTuples(
    gtuples = GTuples(seqnames = mtsv[, chr], 
                      strand = mtsv[, strand], 
                      pos = as.matrix(mtsv[, grep(pattern = '^pos', 
                                                  x = colnames(mtsv)), 
                                           with = FALSE]), 
                      seqinfo = seqinfo), 
    methinfo = methinfo)
  
  # Create the MethPat object
  print("Creating MethPat object...")
  new("MethPat", SummarizedExperiment(assays = assays, rowData = rowData,  
                                      verbose = verbose))
}
