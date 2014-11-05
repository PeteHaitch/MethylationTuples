### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### filter_variants: Filter out variants from MethPat object
###


# TODO: Determine which fields of Bis-SNP VCF are required for filtering.
# TODO: Should filter_variants read in the VCF files or assume the user has 
# done this already and work with the imported objects?

#' Filter out variants from MethPat object.
#' 
#' @param methpat A \code{\link{MethPat}} object.
#' @param variant_files A \code{character} vector with the names of the VCFs 
#' containing the variant calls. Currently only VCFs created by \code{Bis-SNP} 
#' are supported and it is required that there is one VCF per sample and that 
#' the order of \code{variant_files} corresponds to the samples in 
#' \code{methpat}.
#' @param A \code{logical(0)}. If \code{TRUE} rows of the \code{MethPat} object 
#' where all samples have a variant are removed from the object entirely, 
#' otherwise these are retained, albeit with the assay counts set to \code{NA}.
#' @return An updated version of the \code{\link{MethPat}} object, where 
#' variants have had their corresponding assays counts set to \code{NA}. Tuples 
#' where all samples are \code{NA} are retained (unless \code{remove} is 
#' \code{TRUE}).
filter_variants <- function(methpat, variant_files, remove = FALSE, 
                            param = ScanVcfParam(fixed = c("FILTER", ), 
                                                 info = c("CS"))) {
  
  
  # Check param has the necessary parameters
  
  # TODO: readVcf is ignoring the CS field; why?
  
  # Read in the VCFs
  vcfs <- lapply(variant_files, function(file, genome) {
    # TODO: Read in VCF
    vcf <- readVcf(file = file, genome = seqinfo(methpat), )
    # TODO: Convert VCF to minimal GRanges object
    vcf_gr <- 
    # Check that VCF and MethPat are similarly stranded.
    if (.stranded(vcf_gr) != .stranded(methpat)) {
      if (.stranded(methpat)) {
        stop(paste0("'methpat' is stranded, so all 'variant_files' must also ", 
                    "be stranded."))
      } else {
        stop(paste0("'methpat' is unstranded, so all 'variant_files' must ", 
                    "also be unstranded."))
      }
    }
    return(vcf_gr)
  }, methpat = methpat)
  
  # Overlap vcfs and methpat
  
  # NA-ify overlap hits
  
  # Remove NA-rows if remove = TRUE
  
  # Return modified methpat
    
  
  
}
