### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### filterOutVariants: Filter out variants from MethPat object
###

#' Filter out variants from \code{MethPat} object.
#' 
#' @param methpat A \code{\link{MethPat}} object.
#' @param vcfs A named \code{character} vector with the paths of the 
#' \code{VCF}s containing the variant calls. Currently only VCFs created by 
#' \code{Bis-SNP} are supported and it is required that there is one \code{VCF} 
#' per sample and that the order of \code{vcfs} corresponds to the 
#' samples in \code{methpat} (i.e., 
#' \code{identical(names(vcfs), colnames(methpat)} is \code{TRUE}).
#' @param remove A \code{logical(0)}. If \code{TRUE} rows of the \code{MethPat} object 
#' where all samples have a variant are removed from the object entirely, 
#' otherwise these are retained, albeit with the assay counts set to \code{NA}.
#' @param verbose A \code{logical(1)} indicating whether progress 
#' should be reported via calls to \code{\link[base]{message}}.
#' @param bpparam A \code{\link[BiocParallel]{bpparam}} object specifying the 
#' parallelisation strategy, if any. See below for a discussion of 
#' parallelisation options available with \code{readMethtuple}.
#' 
#' @return An updated version of the \code{\link{MethPat}} object, where 
#' variants have had their corresponding assays counts set to \code{NA}. Tuples 
#' where all samples are \code{NA} are retained (unless \code{remove} is 
#' \code{TRUE}).
#' 
#' @section Parallelisation:
#' Parallelisation of \code{filterOutVariants} is partially supported. 
#' \code{VCF}s are read-in and processed in parallel, where appropriate.
#' Parallelisation uses the \pkg{BiocParallel} package. By default this uses a
#' \code{\link[BiocParallel]{MulticoreParam}()} instance or the user's 
#' preferred back-end if they have used \code{\link[BiocParallel]{register}}. 
#' Please consult the \pkg{BiocParallel} documentation for details on 
#' registering a parallel backend and parallelisation support available on 
#' different operating systems.
#' 
#' @export
#' @examples
#' \dontrun{
#' ## TODO
#' }
filterOutVariants <- function(methpat, vcfs, remove = FALSE, 
                                verbose = getOption("verbose"),
                                bpparam = bpparam()) {
  # Check names of vcfs match sample names
  if (!identical(colnames(methpat), names(vcfs))) {
    stop("Names of 'vcfs' must be identical to 'colnames(methpat)'.")
  }
    
  # Read and process the VCFs
  if (verbose) {
    message("Reading and processing VCFs ...")
  }
  vcfs <- bplapply(vcfs, function(vf, methpat) {
    param <- VariantAnnotation::ScanVcfParam(fixed = c("FILTER"), 
                                             info = c("CS"), geno = NA)
    # Read in VCF
    vcf <- VariantAnnotation::readVcf(file = vf, genome = NA_character_, 
                                      param = param)
    
        # Convert VCF to minimal GRanges object
    strand <- VariantAnnotation::info(vcf)$CS
    strand[is.na(strand)] <- '*'
    GRanges(seqnames(vcf), 
            ranges = unname(ranges(vcf)),
            strand = Rle(strand), 
            seqinfo = seqinfo(vcf))[rowRanges(vcf)$FILTER == "PASS"]
  }, methpat = methpat)
  
  # Split m-tuples into list of GRanges storing 1-tuples
  if (verbose) {
    message(paste0("Splitting ", size(methpat), "-tuples ..."))
  }
  mp_grs <- lapply(seq_len(size(methpat)), function(i, seqnames, tuples, 
                                                    strand, seqinfo) {
    GRanges(seqnames = seqnames, 
            ranges = IRanges(tuples[, i], width = 1L), 
            strand = strand, 
            seqinfo = seqinfo)
  }, seqnames = seqnames(methpat), tuples = tuples(methpat), 
  strand = strand(methpat), seqinfo = seqinfo(methpat))
  
  # Identify methylation tuples that overlap
  if (verbose) {
    message(paste0("Identifying ", size(methpat), "-tuples overlapping ", 
                   "variants ..."))
  }
  to_remove <- bplapply(vcfs, function(vcf, mp_grs) {
    # TODO: Investigate creating a GIntervalTree from vcf, which might speed 
    # up the repeated overlap queries against it.
    tr <- lapply(mp_grs, function(mp_gr, vcf) {
      overlapsAny(mp_gr, vcf)
    }, vcf = vcf)
    which(Reduce("|", tr))
  }, mp_grs = mp_grs)
  if (verbose) {
    message(paste0("Replacing counts at ", size(methpat), 
                   "-tuples containing variants with NAs ..."))
    message("Number of loci filtered out per sample:")
    n_remove <- sapply(to_remove, length)
    names(n_remove) <- colnames(methpat)
    print(n_remove)
  }
  # Convert list to vector
  to_remove <- unlist(mapply(function(to_remove, col, nrow) {
    # as.numeric to prevent integer overflow
    to_remove + as.numeric(col) * nrow
  }, to_remove, seq_along(to_remove) - 1L, nrow(methpat)), use.names = FALSE)
  if (!identical(to_remove, numeric(0))) {
    # Can't run in parallel due to usual long-vector problems
    # TODO: Figure out whether I should be using assays(methpat) or 
    # assays(methpat, withDimnames = FALSE)
    assays(methpat) <- SimpleList(lapply(assays(methpat), 
                                         function(assay, to_remove) {
                                           assay[to_remove] <- NA
                                           assay
                                         }, to_remove = to_remove))
  }
  
  # Return modified MethPat object.
  # If remove = TRUE, then remove NA-rows (only need check one assay)
  if (remove) {
    if (verbose) {
      message("Removing NA rows ...")
    }
    return(methpat[rowSums(is.na(assay(methpat, 1))) != ncol(methpat)])
  } else {
    methpat 
  }
}
