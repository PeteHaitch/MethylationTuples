#' Tools for analysing methylation patterns at genomic tuples.
#'
#' \pkg{MethylationTuples} provides tools for analysing, managing and 
#' visualising methylation patterns at genomic tuples. These include analyses 
#' of co-methylation and epipolymorphism.
#'
#' Please refer to the vignettes to see how to use the \pkg{MethylationTuples}
#' package.
#'
#' @docType package
#' @name MethylationTuples-package
#' @useDynLib MethylationTuples, .registration = TRUE
#' @import GenomicTuples
#' @import Rcpp
#' @import BiocGenerics
#' @import IRanges
#' @import data.table
#' @importFrom Biobase validMsg
#' @importFrom BiocParallel bplapply
#' @importFrom R.utils gunzip
NULL
