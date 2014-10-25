### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Internal functions for getting indices to create pairs of methylation loci 
### These are used by betaCor()
###

# TODO: Move as much as possible to C++ and use iterators rather than explicit
# index vectors. This is ridiculously slow for what it does!
#' @export
#' @keywords internal
.makeNeighbourPairsIdx <- function(methpat_rd, mtuples) {
  # x and y index the first and second methlyation loci in each pair, 
  # respectively. The index is with respect to mtuples.
  n <- length(mtuples)
  x <- seq.int(from = 1, to = n - 1, by = 1)
  y <- seq.int(from = 2, to = n, by = 1)
  # Construct indices of the "valid pairs" of methylation loci.
  # A pair is valid if both methylation loci are on the same seqlevel and 
  # strand and the pair.
  # NOTE: It is _much_ faster to extract-then-subset than subset-then-extract,
  # e.g. start(mtuples)[y] is much faster than start(mtuples[y]).
  same_chrom <- seqnames(mtuples)[x] == seqnames(mtuples)[y]
  same_strand <- strand(mtuples[x]) == strand(mtuples[y])
  
  pairs_idx <- which(same_chrom & same_strand)
  pairs <- GRanges(
    seqnames = seqnames(mtuples)[pairs_idx], 
    ranges = IRanges(start = start(mtuples)[x[pairs_idx]],
                     end = start(mtuples)[y[pairs_idx]]),
    strand = strand(mtuples)[x[pairs_idx]])
  
  # Find which valid pairs can be constructed from methpat_rd.
  # Do this by matching methpat_rd to start and end of pairs.
  # Then, remove those loci where only the start or end matches.
  # xx and yy index the first and second methlyation loci in each pair, 
  # respectively, with respect to methpat_rd.
  s <- findOverlaps(query = methpat_rd, subject = pairs, type = 'start')
  e <- findOverlaps(query = methpat_rd, subject = pairs, type = 'end') 
  xx <- queryHits(s)[!is.na(match(subjectHits(s), subjectHits(e)))]
  yy <- queryHits(e)[!is.na(match(subjectHits(e), subjectHits(s)))]
  
  list(x = xx, y = yy)
}

# TODO: Move as much as possible to C++ and use iterators rather than explicit
# index vectors.
# TODO: Cpp_MethylationTuples_makeAllPairs uses obscene amounts of memory 
# because it will make many pairs that are later thrown away, since they are not
# on the same chromosome or have the same strand. Do this smarter.
#' @export
#' @keywords internal
.makeAllPairsIdx <- function(methpat_rd, max_ipd) {
  # p is a list of candidate pairs.
  # p$xx (resp. p$yy) indexes the first (resp. last) element of each pair. 
  # where the index is with respect to methpat_rd.
  p <- .Call(Cpp_MethylationTuples_makeAllPairs, 
             start(methpat_rd),
             max_ipd
  )
  
  # Construct indices of the "valid pairs" of methylation loci.
  # A pair is valid if both methylation loci are on the same seqlevel and 
  # strand and the pair.
  # NOTE: It is _much_ faster to extract-then-subset than subset-then-extract,
  # e.g. start(mtuples)[y] is much faster than start(mtuples[y]).
  same_seqlevel <- seqnames(methpat_rd)[p$x] == seqnames(methpat_rd)[p$y]
  same_strand <- strand(methpat_rd[p$x]) == strand(methpat_rd[p$y])
  # z_vp_idx is the index of valid pairs in (p$x, p$y).
  pairs_idx <- which(same_seqlevel & same_strand)
  
  list(x = p$x[pairs_idx], y = p$y[pairs_idx])   
}