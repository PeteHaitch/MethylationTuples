# TODO: Should this be a function or a method? How to choose when to use a 
# function vs. writing a generic and method?
#' Collapse a \code{\link{MethPat}} object by strand.
#' 
#' A \code{\link{MethPat}} object may contain separate counts for the positive 
#' (\code{+}) and negative (\code{-}) strands. In some instances it may be 
#' desirable to collapse these observations across strands. This function will 
#' combine the loci and sum their associated counts. The co-ordinates of the 
#' collapsed loci will have the \code{\link[GenomicTuples]{strand}} set to 
#' \code{*}, although positions will be with respect to the forward strand. For 
#' example, the CpG with co-ordinates \code{chr1:10} on the \code{+} strand and 
#' \code{chr1:11} on the \code{-} strand will be collapsed to \code{chr1:10} 
#' with the strand set to \code{*}.
#' 
#' \strong{NOTE}: Collapsing by strand is only possible if the 
#' \code{\link{methtype}} of the object is \code{CG}, since neither \code{CHG} 
#' nor \code{CHH} is a strand-symmetric motif.
#' 
#' \strong{NOTE}: Any \code{mcols} of the \code{methpat} object will be dropped.
#' 
#' @param methpat A \code{\link{MethPat}} object.
#' 
#' @return An updated version of the \code{\link{MethPat}} object, where the 
#' loci and their associated counts have been combined across strands.
#' @export
collapseStrand <- function(methpat) {
  if (!is(methpat, "MethPat")) {
    stop("Only defined for MethPat objects.")
  }
  if (!identical(methtype(methpat), "CG")) {
    stop("Can only collapse by strand if 'methtype' is 'CG'.")
  }
  if (!.stranded(methpat)) {
    stop(paste0("'methpat' must be stranded, that is, the 'strand' of all ", 
                "loci must be '+' or '-' and not '*'."))
  }
  if (isTRUE(any(duplicated(methpat)))) {
    stop("'methpat' must not contain any duplicate genomic tuples.")
  }
  # OB_STRAND_OFFSET is the number of base pairs a locus on the '-' (OB) strand 
  # needs to be shifted in order to have the same position as its counterpart 
  # on the '+' (OT) strand. 
  if (identical(methtype(methpat), "CG")) {
    OB_STRAND_OFFSET <- -1L
  }
  
  # Find overlaps between positive and negative strand loci and create indices.
  # TODO: Use gtuples(methpat, use.mcols = FALSE), rather than rowData(methpat), 
  # once implemented.
  rd <- rowData(methpat)

  ol <- suppressWarnings(findOverlaps(rd, shift(rd, OB_STRAND_OFFSET), 
                                      type = "equal", ignore.strand = TRUE))
  
  both <- strand(rd)[queryHits(ol)] == "+" &
    strand(rd)[subjectHits(ol)] == "-"
  plus_both <- queryHits(ol)[as.logical(both)]
  neg_both <- subjectHits(ol)[as.logical(both)]
  plus_only <- setdiff(which(strand(rd) == "+"), plus_both)
  neg_only <- setdiff(which(strand(rd) == "-"), neg_both)
  row_idx <- c(plus_both, plus_only, neg_only)
  
  # Construct new GTuples
  new_rd <- unstrand(rd[row_idx])
  mcols(new_rd) <- NULL
  new_assays <- endoapply(assays(methpat), 
                          function(assay, plus_both, neg_both, plus_only, 
                                   neg_only) {
                            rbind(assay[plus_both, , drop = FALSE] + 
                                    assay[neg_both, , drop = FALSE], 
                                  assay[plus_only, , drop = FALSE],
                                  assay[neg_only, , drop = FALSE])
                            
                          }, plus_both = plus_both, neg_both = neg_both, 
                          plus_only = plus_only, neg_only = neg_only)
  
  MethPat(assays = new_assays, rowData = new_rd, colData = colData(methpat), 
          exptData = exptData(methpat))
}