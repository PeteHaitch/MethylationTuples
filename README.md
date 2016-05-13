# MethylationTuples

[![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/0.1.0/wip.svg)](http://www.repostatus.org/#wip) [![](http://www.r-pkg.org/badges/version/reprex)](http://www.r-pkg.org/pkg/reprex)
[![Linux Build Status](https://travis-ci.org/PeteHaitch/MethylationTuples.svg?branch=DSArray)](https://travis-ci.org/PeteHaitch/MethylationTuples)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/github/PeteHaitch/MethylationTuples?svg=true)](https://ci.appveyor.com/project/PeteHaitch/MethylationTuples)
[![Coverage Status](https://img.shields.io/codecov/c/github/PeteHaitch/MethylationTuples/DSArray.svg)](https://codecov.io/github/PeteHaitch/MethylationTuples?branch=DSArray)

---

_MethylationTuples_ provides tools for analysing, managing and visualising 
methylation patterns at genomic tuples. These include analyses of 
co-methylation and epipolymorphisms.

_MethylationTuples_ is in development and can only be installed using the 
development version of Bioconductor. Please first read 
[these instructions on installing the development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). 

---

## Installation

```r
# Install devtools if not already installed
if (suppressWarnings(!require(devtools))) {
  install.packages("devtools")  
}
# Install development version of MethylationTuples
devtools::install_github("PeteHaitch/MethylationTuples@DSArray")
```

## License

Artistic-2.0