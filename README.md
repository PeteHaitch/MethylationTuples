MethylationTuples
================================================================================
[![Build Status](https://travis-ci.org/PeteHaitch/MethylationTuples.png?branch=master)](https://travis-ci.org/PeteHaitch/MethylationTuples)

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
