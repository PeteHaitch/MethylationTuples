### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### filterOutVariants: Filter out variants from MethPat object
###

# TODO: Clarify NAMESPACE issue with VariantAnnotation, e.g., readVcf doesn't 
# work but VariantAnnotation::readVcf does but I'm not sure this is the correct 
# way of doing things (see http://r-pkgs.had.co.nz/namespace.html)

#' Filter out variants from \code{MethPat} object.
#' 
#' @param methpat A \code{\link{MethPat}} object.
#' @param variant_files A named \code{character} vector with the paths of the 
#' \code{VCF}s containing the variant calls. Currently only VCFs created by 
#' \code{Bis-SNP} are supported and it is required that there is one \code{VCF} 
#' per sample and that the order of \code{variant_files} corresponds to the 
#' samples in \code{methpat} (i.e., 
#' \code{identical(names(variant_files), colnames(methpat)} is \code{TRUE}).
#' @param remove A \code{logical(0)}. If \code{TRUE} rows of the \code{MethPat} object 
#' where all samples have a variant are removed from the object entirely, 
#' otherwise these are retained, albeit with the assay counts set to \code{NA}.
#' @param verbose A \code{\link{logical(1)}} indicating whether progress 
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
#' \code{\link[BiocParallel]{MulticoreParam()}} instance or the user’s 
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
filterOutVariants <- function(methpat, variant_files, remove = FALSE, 
                                verbose = getOption("verbose"),
                                bpparam = bpparam()) {
  # Check names variant_files match sample names
  if (!identical(colnames(methpat), names(variant_files))) {
    stop("Names of 'variant_files' must be identical to 'colnames(methpat)'.")
  }
    
  # Read and process the VCFs
  if (verbose) {
    message("Reading and processing VCFs ...")
  }
  vcfs <- bplapply(variant_files, function(vf, methpat) {
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
            seqinfo = seqinfo(vcf))[rowData(vcf)$FILTER == "PASS"]
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
  to_remove <- do.call("cbind", bplapply(vcfs, function(vcf, mp_grs) {
    # TODO: Investigate creating a GIntervalTree from vcf, which might speed 
    # up the repeated overlap queries against it.
    tr <- lapply(mp_grs, function(mp_gr, vcf) {
      overlapsAny(mp_gr, vcf)
    }, vcf = vcf)
    Reduce("|", tr)
  }, mp_grs = mp_grs))

  # NA-ify overlap hits
  if (verbose) {
    message(paste0("Replacing counts at ", size(methpat), 
                   "-tuples containing variants with NAs ..."))
  }
  assays <- SimpleList(lapply(assays(methpat), function(assay, to_remove) {
    assay[to_remove] <- NA
    assay
  }, to_remove = to_remove))
  assays(methpat) <- assays
  
  # Return modified MethPat object.
  # If remove = TRUE, then remove NA-rows (only need check one assay)
  if (remove) {
    if (verbose) {
      message("Removing NA rows ...")
    }
    return(methpat[rowSums(is.na(assays[[1]])) != ncol(methpat)])
  } else {
    methpat 
  }
}
