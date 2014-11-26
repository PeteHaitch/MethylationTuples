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
  # TODO: names(methpat@assays$field("data")) is 1000x faster than 
  # names(assays(methpat)).
  if (!identical(names(methpat@assays$field("data")), 
                 .makeMethPatNames(size(methpat)))) {
    stop("'methpat' cannot contain non-standard assays.")
  }
  # OB_STRAND_OFFSET is the number of base pairs a locus on the '-' (OB) strand 
  # needs to be shifted in order to have the same position as its counterpart 
  # on the '+' (OT) strand. 
  if (identical(methtype(methpat), "CG")) {
    OB_STRAND_OFFSET <- -1L
  }

  # Find equal overlaps between positive and negative strand loci and create 
  # indices. Not using 
  # findOverlaps(rd, shift(rd, -1L), type = "equal", ignore.strand = TRUE) 
  # because of the crazy memory usage required.
  # This code assumes that there are no duplicate tuples
  # TODO: Use data.table::foverlaps with type = 'equal' once it is implemented.
  rd <- rowData(methpat)
  rd_order <- order(rd)
  # TODO: Shouldn't have to do as.vector(), Rle should just work, I think.
  rd_plus_order <- rd_order[as.vector(strand(rd) == "+")]
  rd_neg_order <- rd_order[as.vector(strand(rd) == "-")]
  rd_plus <- rd[strand(rd) == "+"]
  rd_neg <- rd[strand(rd) == "-"]
  rd_plus_dt <- data.table(seqnames = as.vector(seqnames(rd_plus)), 
                           tuples(rd_plus), plus = rd_plus_order)
  setkeyv(rd_plus_dt, colnames(rd_plus_dt))
  rd_neg_dt <- data.table(seqnames = as.vector(seqnames(rd_neg)), 
                       tuples(rd_neg) + OB_STRAND_OFFSET, neg = rd_neg_order)
  setkeyv(rd_neg_dt, colnames(rd_neg_dt))
  ol <- merge(rd_plus_dt, rd_neg_dt, all = TRUE)
  plus_both <- na.omit(ol[!is.na(neg), plus])
  neg_both <- na.omit(ol[!is.na(plus), neg])
  plus_only <- ol[is.na(neg), plus]
  neg_only <- ol[is.na(plus), neg]
  row_idx <- c(plus_both, plus_only, neg_only)
  
  # Construct new GTuples
  new_rd <- c(unstrand(rd[plus_both]), 
              unstrand(rd[plus_only]),
              shift(unstrand(rd[neg_only]), OB_STRAND_OFFSET))
  mcols(new_rd) <- NULL
  # Construct new assays
  new_assays <- endoapply(assays(methpat), 
                          function(assay, plus_both, neg_both, plus_only, 
                                   neg_only) {
                            # Need special handling of samples where there is 
                            # NA for one strand but not the other. Otherwise 
                            # end up with things like 0 + NA = NA and 
                            # 10 + NA = NA, when I want 0 + NA = 0 and 
                            # 10 + NA = 10.
                            assay_plus <- assay[plus_both, , drop = FALSE]
                            assay_neg <- assay[neg_both, , drop = FALSE]
                            assay_plus[is.na(assay_plus) & 
                                         !is.na(assay_neg)] <- 0L
                            assay_neg[is.na(assay_neg) & 
                                        !is.na(assay_plus)] <- 0L
                            rbind(assay_plus + assay_neg, 
                                  assay[plus_only, , drop = FALSE],
                                  assay[neg_only, , drop = FALSE])
                          }, plus_both = plus_both, neg_both = neg_both, 
                          plus_only = plus_only, neg_only = neg_only)
  # Construct new MethPat object
  MethPat(assays = new_assays, rowData = new_rd, colData = colData(methpat), 
          exptData = exptData(methpat))
}