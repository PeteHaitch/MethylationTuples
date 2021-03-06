% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/findMTuples.R
\name{findMTuples}
\alias{findMTuples}
\title{Find tuples of methylation loci of a given size in a reference genome.}
\usage{
findMTuples(bsgenome, methinfo, size, exclude = NULL)
}
\arguments{
\item{bsgenome}{A \code{\link[BSgenome]{BSgenome}} object for the reference
genome of interest.}

\item{methinfo}{A \code{\link{MethInfo}} object containing the type of
methylation loci for which to search.}

\item{size}{An \code{integer} specifying the size of tuples for which to
search.}

\item{exclude}{A character vector of \code{seqnames} to be filtered out from
the \code{bsgenome}.}
}
\value{
An \code{\link{MTuples}} object containing the tuples.
}
\description{
Find tuples of methylation loci of a given size in a reference genome.
}
\examples{
\dontrun{
# Find all 2-tuples of CG in the C. elegans reference genome
library(BSgenome.Celegans.UCSC.ce2)
findMTuples(bsgenome = Celegans, methinfo = MethInfo('CG'), size = 2)
}
}

