# R package: MethylationTuples

`MethylationTuples` provides tools for analysing, managing and visualising 
methylation patterns at genomic tuples. These include analyses of 
co-methylation and epipolymorphisms.

`MethylationTuples` is in development and can only be installed using the 
development version of Bioconductor. Please first read 
[these instructions on installing the development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). 

```R
# Install devtools if not already installed
if (suppressWarnings(!require(devtools))) {
  install.packages('devtools')  
}
# Install development version of MethylationTuples
devtools::install_github("PeteHaitch/MethylationTuples")
```

## R CMD check status
Travis CI: <a href="https://travis-ci.org/PeteHaitch/MethylationTuples"><img src="https://travis-ci.org/PeteHaitch/MethylationTuples.svg?branch=master" alt="Build status"></a>

## Test coverage status
coveralls.io: [![Coverage Status](https://coveralls.io/repos/PeteHaitch/MethylationTuples/badge.svg?branch=master)](https://coveralls.io/r/PeteHaitch/MethylationTuples?branch=master)

