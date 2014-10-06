MethylationTuples
================================================================================
[![Build Status](https://travis-ci.org/PeteHaitch/MethylationTuples.png?branch=master)](https://travis-ci.org/PeteHaitch/MethylationTuples)

`MethylationTuples` provides tools for analysing, managing and visualising 
methylation patterns at genomic tuples. These include analyses of 
co-methylation and epipolymorphisms.

__This package is in early development and requires the use of the development version of Bioconductor ([instructions here](http://bioconductor.org/developers/how-to/useDevel/)).__ It can be installed (provided it is passing the Travis CI build) with:

```R
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite(c('devtools', 'Rcpp', 'GenomicRanges','testthat', 'knitr', 'BiocStyle'))
devtools::install_github("PeteHaitch/GenomicTuples")
devtools::install_github("PeteHaitch/MethylationTuples")
```

Alternatively, you may like to try using the [`packrat` private library](http://rstudio.github.io/packrat/) that is included in this git repository, which includes the necessary packages for using the `MethylationTuples` package.
